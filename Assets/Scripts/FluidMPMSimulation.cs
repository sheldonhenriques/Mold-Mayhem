using UnityEngine;
using System.Collections.Generic;

// Fluid simulation applying MPM concepts from the clay simulation
public class FluidMPMSimulation : MonoBehaviour
{
    [Header("Fluid Properties")]
    [Range(0.1f, 2f)] public float viscosity = 0.3f;
    [Range(0.5f, 3f)] public float density = 1.0f;
    [Range(0.1f, 2f)] public float surfaceTension = 0.8f;
    [Range(0.1f, 1f)] public float restDensity = 0.9f;
    [Range(0.01f, 0.1f)] public float pressureStiffness = 0.05f;
    
    [Header("Simulation Settings")]
    [Range(0.001f, 0.02f)] public float timeStep = 0.005f;
    [Range(2, 8)] public int substeps = 6;
    [Range(0.05f, 0.3f)] public float particleRadius = 0.1f;
    [Range(16, 128)] public int gridResolution = 48;
    
    [Header("Environmental")]
    public Vector3 gravity = new Vector3(0, -9.81f, 0);
    public Vector3 simulationBounds = new Vector3(5, 5, 5);
    
    private List<FluidParticle> particles = new List<FluidParticle>();
    private FluidGrid grid;
    private MeshFilter meshFilter;
    private MeshRenderer meshRenderer;
    
    public class FluidParticle
    {
        public Vector3 position;
        public Vector3 velocity;
        public Vector3 force;
        public float mass;
        public float density;
        public float temperature; // For steam/heating effects
        public Vector3 gradientPressure;
        public float pressure;
        
        public FluidParticle(Vector3 pos, float m = 1f)
        {
            position = pos;
            velocity = Vector3.zero;
            force = Vector3.zero;
            mass = m;
            density = 0f;
            temperature = 20f; // Room temperature
            gradientPressure = Vector3.zero;
            pressure = 0f;
        }
    }
    
    public class FluidGrid
    {
        public int resolution;
        public Vector3 bounds;
        public float cellSize;
        public FluidCell[,,] cells;
        
        public FluidGrid(int res, Vector3 bounds)
        {
            resolution = res;
            this.bounds = bounds;
            cellSize = bounds.x / resolution;
            cells = new FluidCell[resolution, resolution, resolution];
            
            for (int x = 0; x < resolution; x++)
                for (int y = 0; y < resolution; y++)
                    for (int z = 0; z < resolution; z++)
                        cells[x, y, z] = new FluidCell();
        }
        
        public Vector3Int WorldToGrid(Vector3 worldPos)
        {
            Vector3 normalizedPos = (worldPos + bounds * 0.5f) / cellSize;
            return new Vector3Int(
                Mathf.Clamp(Mathf.FloorToInt(normalizedPos.x), 0, resolution - 1),
                Mathf.Clamp(Mathf.FloorToInt(normalizedPos.y), 0, resolution - 1),
                Mathf.Clamp(Mathf.FloorToInt(normalizedPos.z), 0, resolution - 1)
            );
        }
    }
    
    public class FluidCell
    {
        public Vector3 velocity;
        public float mass;
        public Vector3 force;
        public float density;
        public float pressure;
        
        public void Reset()
        {
            velocity = Vector3.zero;
            mass = 0f;
            force = Vector3.zero;
            density = 0f;
            pressure = 0f;
        }
    }
    
    void Start()
    {
        InitializeComponents();
        InitializeGrid();
        InitializeFluidParticles();
    }
    
    void InitializeComponents()
    {
        meshFilter = GetComponent<MeshFilter>();
        if (meshFilter == null) meshFilter = gameObject.AddComponent<MeshFilter>();
        
        meshRenderer = GetComponent<MeshRenderer>();
        if (meshRenderer == null) meshRenderer = gameObject.AddComponent<MeshRenderer>();
        
        // Create water-like material
        meshRenderer.material = new Material(Shader.Find("Standard"));
        meshRenderer.material.color = new Color(0.2f, 0.6f, 0.9f, 0.7f);
        meshRenderer.material.SetFloat("_Mode", 3); // Transparent mode
        meshRenderer.material.SetInt("_SrcBlend", (int)UnityEngine.Rendering.BlendMode.SrcAlpha);
        meshRenderer.material.SetInt("_DstBlend", (int)UnityEngine.Rendering.BlendMode.OneMinusSrcAlpha);
        meshRenderer.material.SetInt("_ZWrite", 0);
        meshRenderer.material.DisableKeyword("_ALPHATEST_ON");
        meshRenderer.material.EnableKeyword("_ALPHABLEND_ON");
        meshRenderer.material.DisableKeyword("_ALPHAPREMULTIPLY_ON");
        meshRenderer.material.renderQueue = 3000;
    }
    
    void InitializeGrid()
    {
        grid = new FluidGrid(gridResolution, simulationBounds);
    }
    
    void InitializeFluidParticles()
    {
        // Create a water cube
        Vector3 center = Vector3.zero;
        int particlesPerAxis = 15;
        float spacing = 0.12f;
        
        for (int x = 0; x < particlesPerAxis; x++)
        {
            for (int y = 0; y < particlesPerAxis; y++)
            {
                for (int z = 0; z < particlesPerAxis; z++)
                {
                    Vector3 offset = new Vector3(
                        (x - particlesPerAxis * 0.5f) * spacing,
                        (y - particlesPerAxis * 0.5f) * spacing + 1f, // Start above ground
                        (z - particlesPerAxis * 0.5f) * spacing
                    );
                    
                    Vector3 pos = center + offset;
                    particles.Add(new FluidParticle(pos, density));
                }
            }
        }
        
        Debug.Log($"Created {particles.Count} fluid particles");
    }
    
    void Update()
    {
        HandleInput();
        
        // Multi-substep integration for stability
        float dt = timeStep / substeps;
        for (int i = 0; i < substeps; i++)
        {
            SimulateFluidStep(dt);
        }
        
        UpdateFluidMesh();
    }
    
    void HandleInput()
    {
        if (Input.GetMouseButton(0))
        {
            Camera cam = Camera.main;
            if (cam != null)
            {
                Ray ray = cam.ScreenPointToRay(Input.mousePosition);
                Vector3 worldPoint = ray.origin + ray.direction * 3f;
                AddFluidForce(worldPoint, Vector3.down * 15f, 0.8f);
            }
        }
    }
    
    void AddFluidForce(Vector3 worldPos, Vector3 force, float radius)
    {
        foreach (var particle in particles)
        {
            float distance = Vector3.Distance(particle.position, worldPos);
            if (distance < radius)
            {
                float influence = 1f - (distance / radius);
                influence = influence * influence;
                particle.force += force * influence;
                particle.temperature += influence * 5f; // Heat from agitation
            }
        }
    }
    
    void SimulateFluidStep(float dt)
    {
        // Reset grid
        for (int x = 0; x < grid.resolution; x++)
            for (int y = 0; y < grid.resolution; y++)
                for (int z = 0; z < grid.resolution; z++)
                    grid.cells[x, y, z].Reset();
        
        // Calculate particle densities
        CalculateFluidDensities();
        
        // Calculate pressure forces
        CalculatePressureForces();
        
        // Apply viscosity forces
        ApplyViscosityForces();
        
        // Particle to Grid transfer
        ParticleToGrid();
        
        // Update grid velocities
        UpdateGridVelocities(dt);
        
        // Grid to Particle transfer
        GridToParticle(dt);
        
        // Update particle positions
        UpdateParticlePositions(dt);
        
        // Apply boundary conditions
        ApplyFluidBoundaryConditions();
        
        // Cool down particles
        foreach (var particle in particles)
        {
            particle.temperature = Mathf.Max(20f, particle.temperature - dt * 10f);
        }
    }
    
    void CalculateFluidDensities()
    {
        float smoothingRadius = particleRadius * 2.5f;
        
        foreach (var particle in particles)
        {
            particle.density = 0f;
            
            foreach (var neighbor in particles)
            {
                float distance = Vector3.Distance(particle.position, neighbor.position);
                if (distance < smoothingRadius)
                {
                    // Poly6 kernel
                    float h = smoothingRadius;
                    float q = distance / h;
                    if (q <= 1f)
                    {
                        float weight = 315f / (64f * Mathf.PI * Mathf.Pow(h, 9)) * 
                                      Mathf.Pow(h * h - distance * distance, 3);
                        particle.density += neighbor.mass * weight;
                    }
                }
            }
            
            // Calculate pressure based on density difference
            particle.pressure = pressureStiffness * (particle.density - restDensity);
        }
    }
    
    void CalculatePressureForces()
    {
        float smoothingRadius = particleRadius * 2.5f;
        
        foreach (var particle in particles)
        {
            Vector3 pressureForce = Vector3.zero;
            
            foreach (var neighbor in particles)
            {
                if (neighbor == particle) continue;
                
                float distance = Vector3.Distance(particle.position, neighbor.position);
                if (distance < smoothingRadius && distance > 0.001f)
                {
                    Vector3 direction = (particle.position - neighbor.position).normalized;
                    
                    // Spiky kernel gradient for pressure
                    float h = smoothingRadius;
                    float q = distance / h;
                    if (q <= 1f)
                    {
                        float weight = -45f / (Mathf.PI * Mathf.Pow(h, 6)) * 
                                      Mathf.Pow(h - distance, 2);
                        
                        float pressureTerm = (particle.pressure + neighbor.pressure) / 
                                           (2f * neighbor.density);
                        pressureForce += neighbor.mass * pressureTerm * weight * direction;
                    }
                }
            }
            
            particle.gradientPressure = pressureForce;
            particle.force += pressureForce;
        }
    }
    
    void ApplyViscosityForces()
    {
        float smoothingRadius = particleRadius * 2.5f;
        
        foreach (var particle in particles)
        {
            Vector3 viscosityForce = Vector3.zero;
            
            foreach (var neighbor in particles)
            {
                if (neighbor == particle) continue;
                
                float distance = Vector3.Distance(particle.position, neighbor.position);
                if (distance < smoothingRadius)
                {
                    // Viscosity laplacian
                    float h = smoothingRadius;
                    float laplacian = 45f / (Mathf.PI * Mathf.Pow(h, 6)) * (h - distance);
                    
                    Vector3 velocityDiff = neighbor.velocity - particle.velocity;
                    viscosityForce += viscosity * neighbor.mass * velocityDiff * 
                                    laplacian / neighbor.density;
                }
            }
            
            particle.force += viscosityForce;
        }
    }
    
    void ParticleToGrid()
    {
        foreach (var particle in particles)
        {
            Vector3Int gridPos = grid.WorldToGrid(particle.position);
            
            // Distribute to nearby cells with trilinear interpolation
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    for (int dz = -1; dz <= 1; dz++)
                    {
                        int gx = gridPos.x + dx;
                        int gy = gridPos.y + dy;
                        int gz = gridPos.z + dz;
                        
                        if (gx >= 0 && gx < grid.resolution &&
                            gy >= 0 && gy < grid.resolution &&
                            gz >= 0 && gz < grid.resolution)
                        {
                            Vector3 cellCenter = new Vector3(gx, gy, gz) * grid.cellSize - 
                                               grid.bounds * 0.5f;
                            float distance = Vector3.Distance(particle.position, cellCenter);
                            float weight = Mathf.Max(0, 1f - distance / grid.cellSize);
                            
                            var cell = grid.cells[gx, gy, gz];
                            cell.mass += particle.mass * weight;
                            cell.velocity += particle.velocity * particle.mass * weight;
                            cell.density += particle.density * weight;
                        }
                    }
                }
            }
        }
        
        // Normalize by mass
        for (int x = 0; x < grid.resolution; x++)
            for (int y = 0; y < grid.resolution; y++)
                for (int z = 0; z < grid.resolution; z++)
                {
                    var cell = grid.cells[x, y, z];
                    if (cell.mass > Mathf.Epsilon)
                    {
                        cell.velocity /= cell.mass;
                        cell.density /= cell.mass;
                    }
                }
    }
    
    void UpdateGridVelocities(float dt)
    {
        for (int x = 0; x < grid.resolution; x++)
        {
            for (int y = 0; y < grid.resolution; y++)
            {
                for (int z = 0; z < grid.resolution; z++)
                {
                    var cell = grid.cells[x, y, z];
                    if (cell.mass > Mathf.Epsilon)
                    {
                        // Apply gravity
                        cell.velocity += gravity * dt;
                        
                        // Apply viscosity damping
                        cell.velocity *= Mathf.Pow(1f - viscosity * 0.5f, dt);
                    }
                }
            }
        }
    }
    
    void GridToParticle(float dt)
    {
        foreach (var particle in particles)
        {
            Vector3Int gridPos = grid.WorldToGrid(particle.position);
            Vector3 newVelocity = Vector3.zero;
            float totalWeight = 0f;
            
            // Interpolate velocity from grid
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    for (int dz = -1; dz <= 1; dz++)
                    {
                        int gx = gridPos.x + dx;
                        int gy = gridPos.y + dy;
                        int gz = gridPos.z + dz;
                        
                        if (gx >= 0 && gx < grid.resolution &&
                            gy >= 0 && gy < grid.resolution &&
                            gz >= 0 && gz < grid.resolution)
                        {
                            Vector3 cellCenter = new Vector3(gx, gy, gz) * grid.cellSize - 
                                               grid.bounds * 0.5f;
                            float distance = Vector3.Distance(particle.position, cellCenter);
                            float weight = Mathf.Max(0, 1f - distance / grid.cellSize);
                            
                            var cell = grid.cells[gx, gy, gz];
                            newVelocity += cell.velocity * weight;
                            totalWeight += weight;
                        }
                    }
                }
            }
            
            if (totalWeight > Mathf.Epsilon)
            {
                newVelocity /= totalWeight;
                
                // FLIP mixing for fluid behavior
                float flipRatio = 0.99f; // Higher FLIP ratio for fluid
                Vector3 picVelocity = newVelocity;
                Vector3 flipVelocity = particle.velocity + (newVelocity - particle.velocity);
                
                particle.velocity = flipRatio * flipVelocity + (1f - flipRatio) * picVelocity;
            }
            
            // Apply external forces
            if (particle.mass > Mathf.Epsilon)
            {
                particle.velocity += particle.force * dt / particle.mass;
            }
            particle.force = Vector3.zero;
        }
    }
    
    void UpdateParticlePositions(float dt)
    {
        foreach (var particle in particles)
        {
            particle.position += particle.velocity * dt;
        }
    }
    
    void ApplyFluidBoundaryConditions()
    {
        Vector3 halfBounds = simulationBounds * 0.5f;
        float restitution = 0.1f; // Low restitution for fluid
        
        foreach (var particle in particles)
        {
            // X boundaries
            if (particle.position.x < -halfBounds.x)
            {
                particle.position.x = -halfBounds.x;
                particle.velocity.x = -particle.velocity.x * restitution;
            }
            if (particle.position.x > halfBounds.x)
            {
                particle.position.x = halfBounds.x;
                particle.velocity.x = -particle.velocity.x * restitution;
            }
            
            // Y boundaries
            if (particle.position.y < -halfBounds.y)
            {
                particle.position.y = -halfBounds.y;
                particle.velocity.y = -particle.velocity.y * restitution;
            }
            if (particle.position.y > halfBounds.y)
            {
                particle.position.y = halfBounds.y;
                particle.velocity.y = -particle.velocity.y * restitution;
            }
            
            // Z boundaries
            if (particle.position.z < -halfBounds.z)
            {
                particle.position.z = -halfBounds.z;
                particle.velocity.z = -particle.velocity.z * restitution;
            }
            if (particle.position.z > halfBounds.z)
            {
                particle.position.z = halfBounds.z;
                particle.velocity.z = -particle.velocity.z * restitution;
            }
        }
    }
    
    void UpdateFluidMesh()
    {
        // Simple metaball-based mesh generation for fluid surface
        List<Vector3> vertices = new List<Vector3>();
        List<int> triangles = new List<int>();
        
        // Create a simplified fluid surface using particle positions
        int gridRes = 16;
        float cellSize = simulationBounds.x / gridRes;
        float[,,] densityField = new float[gridRes, gridRes, gridRes];
        
        // Generate density field
        for (int x = 0; x < gridRes; x++)
        {
            for (int y = 0; y < gridRes; y++)
            {
                for (int z = 0; z < gridRes; z++)
                {
                    Vector3 worldPos = new Vector3(
                        (x - gridRes * 0.5f) * cellSize,
                        (y - gridRes * 0.5f) * cellSize,
                        (z - gridRes * 0.5f) * cellSize
                    );
                    
                    float density = 0f;
                    foreach (var particle in particles)
                    {
                        float distance = Vector3.Distance(worldPos, particle.position);
                        if (distance < particleRadius * 2f)
                        {
                            float influence = 1f - (distance / (particleRadius * 2f));
                            density += influence * influence * influence;
                        }
                    }
                    densityField[x, y, z] = density;
                }
            }
        }
        
        // Generate mesh from density field (simplified marching cubes)
        float isoValue = 0.3f;
        for (int x = 0; x < gridRes - 1; x++)
        {
            for (int y = 0; y < gridRes - 1; y++)
            {
                for (int z = 0; z < gridRes - 1; z++)
                {
                    if (densityField[x, y, z] > isoValue)
                    {
                        Vector3 center = new Vector3(
                            (x + 0.5f - gridRes * 0.5f) * cellSize,
                            (y + 0.5f - gridRes * 0.5f) * cellSize,
                            (z + 0.5f - gridRes * 0.5f) * cellSize
                        );
                        
                        AddFluidCube(vertices, triangles, center, cellSize * 0.9f);
                    }
                }
            }
        }
        
        // Update mesh
        Mesh fluidMesh = new Mesh();
        if (vertices.Count > 0)
        {
            fluidMesh.vertices = vertices.ToArray();
            fluidMesh.triangles = triangles.ToArray();
            fluidMesh.RecalculateNormals();
            fluidMesh.RecalculateBounds();
        }
        
        meshFilter.mesh = fluidMesh;
    }
    
    void AddFluidCube(List<Vector3> vertices, List<int> triangles, Vector3 center, float size)
    {
        int baseIndex = vertices.Count;
        float half = size * 0.5f;
        
        // Add cube vertices
        Vector3[] cubeVerts = {
            center + new Vector3(-half, -half, -half),
            center + new Vector3(half, -half, -half),
            center + new Vector3(half, half, -half),
            center + new Vector3(-half, half, -half),
            center + new Vector3(-half, -half, half),
            center + new Vector3(half, -half, half),
            center + new Vector3(half, half, half),
            center + new Vector3(-half, half, half)
        };
        
        vertices.AddRange(cubeVerts);
        
        // Add triangles for cube faces
        int[] cubeTriangles = {
            0, 2, 1, 0, 3, 2, // front
            1, 6, 5, 1, 2, 6, // right  
            5, 7, 4, 5, 6, 7, // back
            4, 3, 0, 4, 7, 3, // left
            3, 6, 2, 3, 7, 6, // top
            0, 5, 4, 0, 1, 5  // bottom
        };
        
        foreach (int index in cubeTriangles)
        {
            triangles.Add(baseIndex + index);
        }
    }
    
    void OnDrawGizmos()
    {
        // Draw simulation bounds
        Gizmos.color = Color.cyan;
        Gizmos.DrawWireCube(transform.position, simulationBounds);
        
        // Draw particles (limited for performance)
        if (particles != null && particles.Count < 500)
        {
            foreach (var particle in particles)
            {
                // Color based on velocity magnitude
                float speed = particle.velocity.magnitude;
                Gizmos.color = Color.Lerp(Color.blue, Color.white, speed / 10f);
                Gizmos.DrawSphere(particle.position, particleRadius * 0.5f);
            }
        }
    }
}