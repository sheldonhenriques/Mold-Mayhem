using UnityEngine;
using System.Collections.Generic;
using Unity.Collections;
using Unity.Jobs;
using Unity.Mathematics;
using System.Collections.Concurrent;

// VR-Optimized Firmer Clay Simulation with spherical 3D particles
public class VRClaySimulation : MonoBehaviour
{
    [Header("Material Properties (Adjusted for Firmer Clay)")]
    [Range(0.1f, 10f)] public float elasticity = 3f;          // Increased for more elastic resistance
    [Range(0.1f, 2f)] public float plasticity = 0.5f;          // Reduced plasticity for firmer behavior
    [Range(0.01f, 1f)] public float viscosity = 0.35f;         // Slightly lower viscosity
    [Range(1f, 5f)] public float density = 2f;                 // Higher density for more mass
    [Range(1f, 6f)] public float cohesion = 5f;                // Increased cohesion
    [Range(0.5f, 2f)] public float restDensity = 1.2f;         // Adjusted rest density
    [Range(0.5f, 1f)] public float yieldStrength = 0.6f;     // Higher yield strength before plastic deformation
    [Range(1f, 5f)] public float structuralIntegrity = 3f;     // Higher structural integrity for firm support

    [Header("Simulation Settings")]
    [Range(0.001f, 0.02f)] public float timeStep = 0.008f;
    [Range(2, 8)] public int substeps = 4;
    [Range(0.1f, 2f)] public float deformationForce = 1f;
    [Range(0.05f, 0.5f)] public float interactionRadius = 0.15f;
    [Range(0.05f, 0.2f)] public float particleRadius = 0.06f; // Slightly smaller radius for tighter packing

    [Header("Grid Settings")]
    [Range(16, 64)] public int gridResolution = 32;
    public Vector3 simulationBounds = new Vector3(4, 4, 4);

    [Header("VR Interaction")]
    public Transform leftHandController;
    public Transform rightHandController;
    [Range(0.05f, 0.3f)] public float handInteractionRadius = 0.15f;
    [Range(0.1f, 5f)] public float handForceMultiplier = 2f;

    [Header("Visual")]
    public Material clayMaterial;

    [Header("Performance")]
    [Range(1, 10)] public int meshUpdateFrequency = 3;
    [Range(100, 2000)] public int maxParticles = 800;
    [Range(4, 32)] public int sphereResolution = 8; // Lower for VR performance

    // Native arrays for Job System
    private NativeArray<ParticleData> particles;
    private NativeArray<float3> particleForces;
    private NativeArray<GridCell> gridCells;
    private NativeArray<int> spatialHashGrid;
    private NativeArray<int> spatialHashIndices;

    // Double buffering for spatial hash to avoid race conditions
    private NativeArray<int> spatialHashGridBuffer;
    private NativeArray<int> spatialHashIndicesBuffer;

    // Spatial hashing for neighbor finding
    private int spatialHashSize;
    private float cellSize;

    // Mesh generation for VR
    private Mesh clayMesh;
    private MeshFilter meshFilter;
    private MeshRenderer meshRenderer;
    private MeshCollider meshCollider;
    private int frameCount;

    // Sphere mesh template for instancing
    private Mesh sphereTemplate;
    private Vector3[] sphereVertices;
    private int[] sphereTriangles;
    private Vector3[] sphereNormals;

    // Cached transforms for performance
    private Transform cachedTransform;

    // Job handles for proper dependency management
    private JobHandle simulationJobHandle;
    private JobHandle spatialHashJobHandle;
    private bool useBufferA = true; // For double buffering

    // Copy of particle data for safe gizmo drawing
    private Vector3[] particlePositionsCache;
    private float[] particleTemperaturesCache;
    private bool cacheValid = false;

    // VR hand tracking
    private Vector3[] previousHandPositions = new Vector3[2];
    private bool[] handWasTracked = new bool[2];

    [System.Serializable]
    public struct ParticleData
    {
        public float3 position;
        public float3 velocity;
        public float3 restPosition;
        public float mass;
        public float temperature;
        public float density;
        public float plasticDeformation;
        public float3 plasticStrain;
        public int neighborCount;
        public float structuralStress;
    }

    public struct GridCell
    {
        public float3 velocity;
        public float mass;
        public float density;
    }

    // Job for updating spatial hash
    [Unity.Burst.BurstCompile]
    public struct SpatialHashJob : IJob
    {
        [ReadOnly] public NativeArray<ParticleData> particles;
        [WriteOnly] public NativeArray<int> spatialHashGrid;
        [WriteOnly] public NativeArray<int> spatialHashIndices;
        [ReadOnly] public int spatialHashSize;
        [ReadOnly] public float cellSize;

        public void Execute()
        {
            // Clear spatial hash
            for (int i = 0; i < spatialHashSize; i++)
            {
                spatialHashGrid[i] = -1;
            }

            // Build spatial hash
            for (int i = 0; i < particles.Length; i++)
            {
                var particle = particles[i];
                int hash = GetSpatialHash(particle.position);

                spatialHashIndices[i] = hash;
                if (spatialHashGrid[hash] == -1)
                {
                    spatialHashGrid[hash] = i;
                }
            }
        }

        private int GetSpatialHash(float3 position)
        {
            int3 gridPos = new int3(
                (int)math.floor(position.x / cellSize),
                (int)math.floor(position.y / cellSize),
                (int)math.floor(position.z / cellSize)
            );

            const int p1 = 73856093;
            const int p2 = 19349663;
            const int p3 = 83492791;

            int hash = (gridPos.x * p1 ^ gridPos.y * p2 ^ gridPos.z * p3) % spatialHashSize;
            return hash < 0 ? hash + spatialHashSize : hash;
        }
    }

    // Job for particle simulation using Burst compilation
    [Unity.Burst.BurstCompile]
    public struct ParticleSimulationJob : IJobParallelFor
    {
        [ReadOnly] public NativeArray<ParticleData> particles;
        [ReadOnly] public NativeArray<int> spatialHash;
        [ReadOnly] public NativeArray<int> spatialIndices;
        [WriteOnly] public NativeArray<float3> forces;

        [ReadOnly] public float deltaTime;
        [ReadOnly] public float cohesionStrength;
        [ReadOnly] public float particleRadius;
        [ReadOnly] public float elasticity;
        [ReadOnly] public float yieldStrength;
        [ReadOnly] public float structuralIntegrity;
        [ReadOnly] public int spatialHashSize;
        [ReadOnly] public float cellSize;
        [ReadOnly] public float3 simulationBounds;

        public void Execute(int index)
        {
            var particle = particles[index];
            float3 totalForce = float3.zero;
            float3 cohesiveForce = float3.zero;
            float3 structuralForce = float3.zero;
            int neighborCount = 0;

            // Calculate forces from neighbors
            for (int i = 0; i < particles.Length; i++)
            {
                if (i == index) continue;

                var neighbor = particles[i];
                float distance = math.distance(particle.position, neighbor.position);
                float influenceRadius = particleRadius * 2.5f; // Tighter interaction radius for firmer clay

                if (distance < influenceRadius && distance > 0.001f)
                {
                    float3 direction = math.normalize(neighbor.position - particle.position);
                    neighborCount++;

                    // Strong cohesion force for firmer clay behavior
                    float cohForce = cohesionStrength * (1f - distance / influenceRadius);
                    cohesiveForce += direction * cohForce * 0.7f; // Increase weight on cohesive force

                    // Structural integrity - resist deformation strongly
                    float neighborRestDistance = math.distance(particle.restPosition, neighbor.restPosition);
                    float structuralStrain = math.abs(distance - neighborRestDistance) / math.max(neighborRestDistance, 0.001f);

                    if (structuralStrain > yieldStrength)
                    {
                        // Plastic deformation occurs, but limited for firmer clay
                        float plasticResponse = (structuralStrain - yieldStrength) * 0.05f;
                        structuralForce += direction * plasticResponse;
                    }
                    else
                    {
                        // Elastic response - maintain original structure strongly
                        float elasticResponse = (neighborRestDistance - distance) * structuralIntegrity;
                        structuralForce += direction * elasticResponse;
                    }

                    // Strong pressure force to prevent overlap
                    if (distance < particleRadius * 1.8f) // Tighter pressure threshold
                    {
                        float pressureStrength = (particleRadius * 1.8f - distance) / (particleRadius * 1.8f);
                        totalForce -= direction * pressureStrength * 2f; // Increased pressure repulsion
                    }
                }
            }

            // Apply cohesive and structural forces
            totalForce += cohesiveForce;
            totalForce += structuralForce;

            // Strong elastic force towards rest position (key for clay behavior)
            float3 restDirection = particle.restPosition - particle.position;
            float restDistance = math.length(restDirection);
            if (restDistance > 0.001f)
            {
                float elasticStrength = math.min(restDistance * elasticity, 5f); // Increased cap for elastic force
                totalForce += math.normalize(restDirection) * elasticStrength;
            }

            // Reduced gravity for clay
            totalForce += new float3(0, -2.5f, 0) * particle.mass;

            // Damping based on neighbor density (more neighbors = more damping)
            float dampingFactor = math.max(0.2f, 1f - (neighborCount * 0.1f)); // Slightly less damping to keep shape
            totalForce *= dampingFactor;

            forces[index] = totalForce;
        }
    }

    // Job for updating particle positions
    [Unity.Burst.BurstCompile]
    public struct ParticleIntegrationJob : IJobParallelFor
    {
        public NativeArray<ParticleData> particles;
        [ReadOnly] public NativeArray<float3> forces;
        [ReadOnly] public float deltaTime;
        [ReadOnly] public float viscosity;
        [ReadOnly] public float yieldStrength;
        [ReadOnly] public float3 simulationBounds;

        public void Execute(int index)
        {
            var particle = particles[index];

            // Update velocity with clay-specific damping
            if (particle.mass > 0.001f)
            {
                particle.velocity += forces[index] * deltaTime / particle.mass;
            }

            // Increased viscosity for firmer clay - more resistance to flow
            float clayViscosity = viscosity * (1f + particle.structuralStress * 2f);
            particle.velocity *= math.pow(1f - clayViscosity, deltaTime * 0.8f);

            // Limit maximum velocity to prevent explosive behavior
            float maxVelocity = 3f; // Lower max velocity for firmer clay
            float velocityMagnitude = math.length(particle.velocity);
            if (velocityMagnitude > maxVelocity)
            {
                particle.velocity = math.normalize(particle.velocity) * maxVelocity;
            }

            // Update position
            float3 oldPosition = particle.position;
            particle.position += particle.velocity * deltaTime;

            // Calculate plastic deformation
            float3 displacement = particle.position - particle.restPosition;
            float displacementMagnitude = math.length(displacement);

            if (displacementMagnitude > yieldStrength)
            {
                // Permanent deformation occurs, but limited
                float plasticAmount = (displacementMagnitude - yieldStrength) * 0.15f;
                particle.restPosition += math.normalize(displacement) * plasticAmount;
                particle.plasticDeformation += plasticAmount;
                particle.structuralStress = math.min(particle.structuralStress + plasticAmount * 0.3f, 1f);
            }
            else
            {
                // Gradual stress relief
                particle.structuralStress = math.max(0, particle.structuralStress - deltaTime * 0.1f);
            }

            // Boundary conditions with clay-specific behavior
            float3 halfBounds = simulationBounds * 0.5f;
            float damping = 0.5f; // Increased damping for firmer bounce

            if (particle.position.x < -halfBounds.x)
            {
                particle.position.x = -halfBounds.x;
                particle.velocity.x = -particle.velocity.x * damping;
                particle.restPosition.x = math.max(particle.restPosition.x, particle.position.x + 0.05f);
            }
            else if (particle.position.x > halfBounds.x)
            {
                particle.position.x = halfBounds.x;
                particle.velocity.x = -particle.velocity.x * damping;
                particle.restPosition.x = math.min(particle.restPosition.x, particle.position.x - 0.05f);
            }

            if (particle.position.y < -halfBounds.y)
            {
                particle.position.y = -halfBounds.y + 0.01f;
                particle.velocity.y = math.abs(particle.velocity.y) * damping;
                // Clay sticks to ground - adjust rest position
                particle.restPosition.y = math.max(particle.restPosition.y, particle.position.y + 0.1f);
            }
            else if (particle.position.y > halfBounds.y)
            {
                particle.position.y = halfBounds.y;
                particle.velocity.y = -particle.velocity.y * damping;
            }

            if (particle.position.z < -halfBounds.z)
            {
                particle.position.z = -halfBounds.z;
                particle.velocity.z = -particle.velocity.z * damping;
                particle.restPosition.z = math.max(particle.restPosition.z, particle.position.z + 0.05f);
            }
            else if (particle.position.z > halfBounds.z)
            {
                particle.position.z = halfBounds.z;
                particle.velocity.z = -particle.velocity.z * damping;
                particle.restPosition.z = math.min(particle.restPosition.z, particle.position.z - 0.05f);
            }

            // Cool down temperature gradually
            if (particle.temperature > 0.1f)
            {
                particle.temperature = math.max(0, particle.temperature - deltaTime * 0.5f); // Slower cooling
            }

            particles[index] = particle;
        }
    }

    void Start()
    {
        InitializeComponents();
        InitializeNativeArrays();
        InitializeParticles();
        CreateSphereTemplate();
        cachedTransform = transform;

        // Initialize cache arrays
        particlePositionsCache = new Vector3[maxParticles];
        particleTemperaturesCache = new float[maxParticles];

        // Initialize hand tracking
        previousHandPositions[0] = Vector3.zero;
        previousHandPositions[1] = Vector3.zero;
        handWasTracked[0] = false;
        handWasTracked[1] = false;
    }

    void InitializeComponents()
    {
        if (this == null || gameObject == null) return;

        meshFilter = GetComponent<MeshFilter>();
        if (meshFilter == null) meshFilter = gameObject.AddComponent<MeshFilter>();

        meshRenderer = GetComponent<MeshRenderer>();
        if (meshRenderer == null) meshRenderer = gameObject.AddComponent<MeshRenderer>();

        meshCollider = GetComponent<MeshCollider>();
        if (meshCollider == null) meshCollider = gameObject.AddComponent<MeshCollider>();

        if (clayMaterial != null)
            meshRenderer.material = clayMaterial;
        else
        {
            meshRenderer.material = new Material(Shader.Find("Standard"));
            meshRenderer.material.color = new Color(0.6f, 0.5f, 0.3f);
        }
            
        clayMesh = new Mesh();
        clayMesh.name = "VR Firmer Clay Mesh";
        clayMesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        meshFilter.mesh = clayMesh;
    }

    void CreateSphereTemplate()
    {
        // Create a low-poly sphere template for VR performance
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        List<Vector3> normals = new List<Vector3>();

        // Generate sphere vertices using spherical coordinates
        int rings = sphereResolution;
        int sectors = sphereResolution;

        for (int r = 0; r <= rings; r++)
        {
            float v = (float)r / rings;
            float phi = v * Mathf.PI;

            for (int s = 0; s <= sectors; s++)
            {
                float u = (float)s / sectors;
                float theta = u * 2 * Mathf.PI;

                Vector3 vertex = new Vector3(
                    Mathf.Sin(phi) * Mathf.Cos(theta),
                    Mathf.Cos(phi),
                    Mathf.Sin(phi) * Mathf.Sin(theta)
                );

                vertices.Add(vertex);
                normals.Add(vertex.normalized);
            }
        }

        // Generate triangles
        for (int r = 0; r < rings; r++)
        {
            for (int s = 0; s < sectors; s++)
            {
                int current = r * (sectors + 1) + s;
                int next = current + sectors + 1;

                // First triangle
                triangles.Add(current);
                triangles.Add(next);
                triangles.Add(current + 1);

                // Second triangle
                triangles.Add(current + 1);
                triangles.Add(next);
                triangles.Add(next + 1);
            }
        }

        sphereVertices = vertices.ToArray();
        sphereTriangles = triangles.ToArray();
        sphereNormals = normals.ToArray();
    }

    void InitializeNativeArrays()
    {
        // Calculate spatial hash parameters
        cellSize = particleRadius * 2f;
        spatialHashSize = Mathf.NextPowerOfTwo(maxParticles * 2);

        // Initialize native arrays
        particles = new NativeArray<ParticleData>(maxParticles, Allocator.Persistent);
        particleForces = new NativeArray<float3>(maxParticles, Allocator.Persistent);

        // Double buffering for spatial hash
        spatialHashGrid = new NativeArray<int>(spatialHashSize, Allocator.Persistent);
        spatialHashIndices = new NativeArray<int>(maxParticles, Allocator.Persistent);
        spatialHashGridBuffer = new NativeArray<int>(spatialHashSize, Allocator.Persistent);
        spatialHashIndicesBuffer = new NativeArray<int>(maxParticles, Allocator.Persistent);

        // Initialize grid
        int gridSize = gridResolution * gridResolution * gridResolution;
        gridCells = new NativeArray<GridCell>(gridSize, Allocator.Persistent);
    }

    void InitializeParticles()
    {
        Vector3 center = new Vector3(0, 1f, 0);
        float radius = 0.8f;
        int particlesPerAxis = Mathf.RoundToInt(Mathf.Pow(maxParticles, 1f/3f));
        float spacing = particleRadius * 1.6f; // Tighter packing for firmer mass

        int particleCount = 0;

        for (int x = 0; x < particlesPerAxis && particleCount < maxParticles; x++)
        {
            for (int y = 0; y < particlesPerAxis && particleCount < maxParticles; y++)
            {
                for (int z = 0; z < particlesPerAxis && particleCount < maxParticles; z++)
                {
                    Vector3 offset = new Vector3(
                        (x - particlesPerAxis * 0.5f) * spacing,
                        (y - particlesPerAxis * 0.5f) * spacing,
                        (z - particlesPerAxis * 0.5f) * spacing
                    );
                    
                    Vector3 pos = center + offset;
                    
                    if (offset.magnitude <= radius)
                    {
                        particles[particleCount] = new ParticleData
                        {
                            position = pos,
                            velocity = float3.zero,
                            restPosition = pos,
                            mass = density,
                            temperature = 0f,
                            density = 0f,
                            plasticDeformation = 0f,
                            plasticStrain = float3.zero,
                            neighborCount = 0,
                            structuralStress = 0f
                        };
                        particleCount++;
                    }
                }
            }
        }
        
        Debug.Log($"Created {particleCount} VR clay particles (firmer configuration)");
    }

    void Update()
    {
        // Complete all previous jobs before starting new frame
        simulationJobHandle.Complete();
        
        HandleVRInput();
        
        // Run simulation with proper job dependencies
        float dt = timeStep / substeps;
        for (int i = 0; i < substeps; i++)
        {
            SimulateStepJob(dt);
        }
        
        // Update mesh less frequently for VR performance
        frameCount++;
        if (frameCount % meshUpdateFrequency == 0)
        {
            UpdateMeshVR();
        }
    }

    void LateUpdate()
    {
        // Cache particle data for safe gizmo drawing after all jobs are complete
        if (particles.IsCreated && simulationJobHandle.IsCompleted)
        {
            UpdateParticleCache();
        }
    }

    void UpdateParticleCache()
    {
        // Ensure jobs are complete before accessing particle data
        simulationJobHandle.Complete();

        for (int i = 0; i < particles.Length; i++)
        {
            var particle = particles[i];
            particlePositionsCache[i] = particle.position;
            particleTemperaturesCache[i] = particle.temperature;
        }
        cacheValid = true;
    }

    void HandleVRInput()
    {
        // Handle left hand controller
        if (leftHandController != null)
        {
            HandleHandInteraction(leftHandController, 0);
        }
        
        // Handle right hand controller
        if (rightHandController != null)
        {
            HandleHandInteraction(rightHandController, 1);
        }
        
        // Fallback to mouse input for testing
        if (leftHandController == null && rightHandController == null)
        {
            HandleMouseInput();
        }
    }

    void HandleHandInteraction(Transform handTransform, int handIndex)
    {
        Vector3 currentHandPosition = handTransform.position;
        
        if (handWasTracked[handIndex])
        {
            Vector3 handVelocity = (currentHandPosition - previousHandPositions[handIndex]) / Time.deltaTime;
            float handSpeed = handVelocity.magnitude;
            
            // Only apply force if hand is moving with sufficient speed
            if (handSpeed > 0.1f)
            {
                // Complete jobs before modifying particle data
                simulationJobHandle.Complete();
                
                Vector3 force = handVelocity * handForceMultiplier * handSpeed;
                ApplyExternalForceOptimized(currentHandPosition, force, handInteractionRadius);
                
                // Invalidate cache since we modified particle data
                cacheValid = false;
            }
        }
        
        previousHandPositions[handIndex] = currentHandPosition;
        handWasTracked[handIndex] = true;
    }

    void HandleMouseInput()
    {
        if (Input.GetMouseButton(0) || Input.GetMouseButton(1))
        {
            // Complete jobs before modifying particle data
            simulationJobHandle.Complete();
            
            Camera cam = Camera.main ?? FindObjectOfType<Camera>();
            if (cam != null)
            {
                Ray ray = cam.ScreenPointToRay(Input.mousePosition);
                Vector3 worldPoint = ray.origin + ray.direction * 5f;
                float forceMultiplier = Input.GetMouseButton(0) ? -1f : 1f;
                ApplyExternalForceOptimized(worldPoint, ray.direction * deformationForce * forceMultiplier, interactionRadius);
                
                // Invalidate cache since we modified particle data
                cacheValid = false;
            }
        }
    }

    void ApplyExternalForceOptimized(Vector3 worldPos, Vector3 force, float radius)
    {
        float3 worldPos3 = worldPos;
        float3 force3 = force;
        float radiusSqr = radius * radius;
        
        for (int i = 0; i < particles.Length; i++)
        {
            var particle = particles[i];
            float distanceSqr = math.distancesq(particle.position, worldPos3);
            
            if (distanceSqr < radiusSqr)
            {
                float distance = math.sqrt(distanceSqr);
                float influence = 1f - (distance / radius);
                influence = influence * influence;
                
                particle.velocity += force3 * influence * 0.4f; // Reduced influence for firmer clay
                particle.temperature += influence * 0.3f;
                
                // For strong deformation, allow limited plastic deformation
                if (influence > 0.5f)
                {
                    float3 deformationDirection = math.normalize(force3);
                    particle.restPosition += deformationDirection * influence * 0.03f;
                    particle.plasticDeformation += influence * 0.1f;
                }
                particles[i] = particle;
            }
        }
    }

    void SimulateStepJob(float dt)
    {
        // Complete previous jobs before starting new ones
        simulationJobHandle.Complete();
        spatialHashJobHandle.Complete();
        
        // Simplified approach: Skip spatial hashing for now to avoid dependency issues
        var forceJob = new ParticleSimulationJob
        {
            particles = particles,
            spatialHash = spatialHashGrid, // Not used in simplified version
            spatialIndices = spatialHashIndices, // Not used in simplified version
            forces = particleForces,
            deltaTime = dt,
            cohesionStrength = cohesion,
            particleRadius = particleRadius,
            elasticity = elasticity,
            yieldStrength = yieldStrength,
            structuralIntegrity = structuralIntegrity,
            spatialHashSize = spatialHashSize,
            cellSize = cellSize,
            simulationBounds = simulationBounds
        };
        
        var forceJobHandle = forceJob.Schedule(particles.Length, 32);
        
        // Schedule integration job (depends on force job)
        var integrationJob = new ParticleIntegrationJob
        {
            particles = particles,
            forces = particleForces,
            deltaTime = dt,
            viscosity = viscosity,
            yieldStrength = yieldStrength,
            simulationBounds = simulationBounds
        };
        
        simulationJobHandle = integrationJob.Schedule(particles.Length, 32, forceJobHandle);
        
        // Mark cache as invalid since particles will be modified
        cacheValid = false;
    }

    void UpdateMeshVR()
    {
        simulationJobHandle.Complete();
        
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        List<Vector3> normals = new List<Vector3>();
        
        CreateSphericalParticleMesh(vertices, triangles, normals);
        
        if (vertices.Count > 0)
        {
            clayMesh.Clear();
            clayMesh.vertices = vertices.ToArray();
            clayMesh.triangles = triangles.ToArray();
            clayMesh.normals = normals.ToArray();
            clayMesh.RecalculateBounds();
            
            if (meshCollider != null)
            {
                meshCollider.sharedMesh = null;
                meshCollider.sharedMesh = clayMesh;
            }
        }
    }

    void CreateSphericalParticleMesh(List<Vector3> vertices, List<int> triangles, List<Vector3> normals)
    {
        for (int i = 0; i < particles.Length; i++)
        {
            var particle = particles[i];
            if (!math.all(math.isfinite(particle.position)))
                continue;
            Vector3 particlePos = particle.position;
            float scale = particleRadius;
            
            // Adjust scale based on temperature (heated particles appear larger)
            scale *= 1f + particle.temperature * 0.1f; // Reduced temperature influence for firmer shape
            
            int baseVertexIndex = vertices.Count;
            
            // Add sphere vertices for this particle
            for (int v = 0; v < sphereVertices.Length; v++)
            {
                Vector3 scaledVertex = sphereVertices[v] * scale + particlePos;
                vertices.Add(scaledVertex);
                normals.Add(sphereNormals[v]);
            }
            
            // Add sphere triangles for this particle
            for (int t = 0; t < sphereTriangles.Length; t++)
            {
                triangles.Add(sphereTriangles[t] + baseVertexIndex);
            }
        }
    }
    void OnDestroy()
    {
        // Complete all jobs before disposing
        simulationJobHandle.Complete();

        // Clean up native arrays
        if (particles.IsCreated) particles.Dispose();
        if (particleForces.IsCreated) particleForces.Dispose();
        if (spatialHashGrid.IsCreated) spatialHashGrid.Dispose();
        if (spatialHashIndices.IsCreated) spatialHashIndices.Dispose();
        if (spatialHashGridBuffer.IsCreated) spatialHashGridBuffer.Dispose();
        if (spatialHashIndicesBuffer.IsCreated) spatialHashIndicesBuffer.Dispose();
        if (gridCells.IsCreated) gridCells.Dispose();
        
        if (clayMesh != null)
        {
            DestroyImmediate(clayMesh);
            clayMesh = null;
        }
    }
    
    void OnDrawGizmos()
    {
        // Use cached data to avoid job system conflicts
        if (this == null || gameObject == null) return;
        if (!cacheValid || particlePositionsCache == null) return;
        
        Gizmos.color = Color.yellow;
        if (cachedTransform != null)
            Gizmos.DrawWireCube(cachedTransform.position, simulationBounds);
        
        // Draw only every 10th particle for performance
        for (int i = 0; i < particlePositionsCache.Length && i < maxParticles; i += 10)
        {
            if (i < particleTemperaturesCache.Length)
            {
                Gizmos.color = Color.Lerp(Color.blue, Color.red, particleTemperaturesCache[i]);
                Gizmos.DrawSphere(particlePositionsCache[i], particleRadius * 0.5f);
            }
        }
        
        // Draw VR controller interaction spheres
        if (leftHandController != null)
        {
            Gizmos.color = Color.green;
            Gizmos.DrawWireSphere(leftHandController.position, handInteractionRadius);
        }
        
        if (rightHandController != null)
        {
            Gizmos.color = Color.blue;
            Gizmos.DrawWireSphere(rightHandController.position, handInteractionRadius);
        }
    }
}