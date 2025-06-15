using UnityEngine;
using System.Collections.Generic;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class ClayPhysics : MonoBehaviour
{
    [Header("Clay Properties")]
    [Range(0.1f, 2.0f)]
    public float elasticity = 0.8f;
    
    [Range(0.0f, 1.0f)]
    public float plasticity = 0.3f;
    
    [Range(0.1f, 5.0f)]
    public float viscosity = 1.0f;
    
    [Range(0.01f, 0.5f)]
    public float particleRadius = 0.05f;
    
    [Range(10, 1000)]
    public int particleCount = 200;
    
    [Header("Deformation Settings")]
    public float maxDeformation = 0.5f;
    public float recoverySpeed = 2.0f;
    public LayerMask handLayer = -1;
    
    [Header("VR Hand References")]
    public Transform leftHand;
    public Transform rightHand;
    public float handInteractionRadius = 0.1f;
    
    // Internal variables
    private ClayParticle[] particles;
    private Vector3[] originalVertices;
    private Vector3[] currentVertices;
    private MeshFilter meshFilter;
    private Mesh clayMesh;
    private Bounds clayBounds;
    
    // Performance optimization
    private float lastUpdateTime;
    private const float UPDATE_INTERVAL = 0.016f; // ~60fps
    
    [System.Serializable]
    public struct ClayParticle
    {
        public Vector3 position;
        public Vector3 velocity;
        public Vector3 force;
        public float mass;
        public float pressure;
        public float density;
        public bool isDeformed;
        public float deformationAmount;
        public Vector3 originalPosition;
    }
    
    void Start()
    {
        InitializeClaySystem();
    }
    
    void InitializeClaySystem()
    {
        meshFilter = GetComponent<MeshFilter>();
        clayMesh = Instantiate(meshFilter.sharedMesh);
        meshFilter.mesh = clayMesh;
        
        // Store original vertices
        originalVertices = clayMesh.vertices;
        currentVertices = new Vector3[originalVertices.Length];
        System.Array.Copy(originalVertices, currentVertices, originalVertices.Length);
        
        // Initialize clay bounds
        clayBounds = GetComponent<Renderer>().bounds;
        
        // Create particle system
        InitializeParticles();
        
        Debug.Log($"Clay Physics initialized with {particleCount} particles");
    }
    
    void InitializeParticles()
    {
        particles = new ClayParticle[particleCount];
        
        for (int i = 0; i < particleCount; i++)
        {
            particles[i] = new ClayParticle
            {
                position = GetRandomPositionInBounds(),
                velocity = Vector3.zero,
                force = Vector3.zero,
                mass = 1.0f,
                pressure = 0.0f,
                density = 1.0f,
                isDeformed = false,
                deformationAmount = 0.0f
            };
            
            particles[i].originalPosition = particles[i].position;
        }
    }
    
    Vector3 GetRandomPositionInBounds()
    {
        Vector3 localPos = new Vector3(
            Random.Range(-0.5f, 0.5f),
            Random.Range(-0.5f, 0.5f),
            Random.Range(-0.5f, 0.5f)
        );
        
        return transform.TransformPoint(localPos);
    }
    
    void Update()
    {
        if (Time.time - lastUpdateTime >= UPDATE_INTERVAL)
        {
            UpdateClayPhysics();
            lastUpdateTime = Time.time;
        }
    }
    
    void UpdateClayPhysics()
    {
        // Reset forces
        for (int i = 0; i < particles.Length; i++)
        {
            particles[i].force = Vector3.zero;
        }
        
        // Apply elastoplasticity constraints
        ApplyElastoplasticityConstraints();
        
        // Handle VR hand interactions
        HandleVRHandInteractions();
        
        // Update particle positions
        UpdateParticlePositions();
        
        // Update mesh deformation
        UpdateMeshDeformation();
        
        // Apply recovery
        ApplyRecovery();
    }
    
    void ApplyElastoplasticityConstraints()
    {
        for (int i = 0; i < particles.Length; i++)
        {
            for (int j = i + 1; j < particles.Length; j++)
            {
                Vector3 displacement = particles[j].position - particles[i].position;
                float distance = displacement.magnitude;
                
                if (distance < particleRadius * 2.0f && distance > 0.001f)
                {
                    Vector3 direction = displacement.normalized;
                    float overlap = (particleRadius * 2.0f) - distance;
                    
                    // Elastic force
                    Vector3 elasticForce = direction * overlap * elasticity;
                    
                    // Plastic deformation check
                    float strain = overlap / (particleRadius * 2.0f);
                    if (strain > plasticity)
                    {
                        // Apply plastic deformation
                        particles[i].isDeformed = true;
                        particles[j].isDeformed = true;
                        particles[i].deformationAmount = Mathf.Min(strain, maxDeformation);
                        particles[j].deformationAmount = Mathf.Min(strain, maxDeformation);
                        
                        // Reduce elastic force for plastic behavior
                        elasticForce *= (1.0f - plasticity);
                    }
                    
                    // Apply viscosity damping
                    Vector3 relativeVelocity = particles[j].velocity - particles[i].velocity;
                    Vector3 viscousForce = relativeVelocity * viscosity * -1.0f;
                    
                    // Apply forces
                    particles[i].force -= elasticForce + viscousForce;
                    particles[j].force += elasticForce + viscousForce;
                }
            }
        }
    }
    
    void HandleVRHandInteractions()
    {
        HandleHandInteraction(leftHand);
        HandleHandInteraction(rightHand);
    }
    
    void HandleHandInteraction(Transform hand)
    {
        if (hand == null) return;
        
        Vector3 handWorldPos = hand.position;
        
        for (int i = 0; i < particles.Length; i++)
        {
            float distance = Vector3.Distance(particles[i].position, handWorldPos);
            
            if (distance < handInteractionRadius)
            {
                // Calculate interaction force
                Vector3 direction = (particles[i].position - handWorldPos).normalized;
                float forceMagnitude = (handInteractionRadius - distance) / handInteractionRadius;
                
                // Apply deformation
                particles[i].isDeformed = true;
                particles[i].deformationAmount = Mathf.Max(particles[i].deformationAmount, forceMagnitude * maxDeformation);
                
                // Push particles away from hand
                Vector3 pushForce = direction * forceMagnitude * 10.0f;
                particles[i].force += pushForce;
                
                // Add some hand velocity influence
                if (hand.GetComponent<Rigidbody>() != null)
                {
                    particles[i].velocity += hand.GetComponent<Rigidbody>().linearVelocity * 0.1f;
                }
            }
        }
    }
    
    void UpdateParticlePositions()
    {
        float deltaTime = Time.deltaTime;
        
        for (int i = 0; i < particles.Length; i++)
        {
            // Integrate velocity
            particles[i].velocity += (particles[i].force / particles[i].mass) * deltaTime;
            
            // Apply air resistance
            particles[i].velocity *= 0.98f;
            
            // Integrate position
            particles[i].position += particles[i].velocity * deltaTime;
            
            // Keep particles within bounds
            ConstrainToBounds(ref particles[i]);
        }
    }
    
    void ConstrainToBounds(ref ClayParticle particle)
    {
        Vector3 localPos = transform.InverseTransformPoint(particle.position);
        
        if (localPos.x < -0.6f || localPos.x > 0.6f ||
            localPos.y < -0.6f || localPos.y > 0.6f ||
            localPos.z < -0.6f || localPos.z > 0.6f)
        {
            // Bounce back into bounds
            if (localPos.x < -0.6f) { localPos.x = -0.6f; particle.velocity.x *= -0.5f; }
            if (localPos.x > 0.6f) { localPos.x = 0.6f; particle.velocity.x *= -0.5f; }
            if (localPos.y < -0.6f) { localPos.y = -0.6f; particle.velocity.y *= -0.5f; }
            if (localPos.y > 0.6f) { localPos.y = 0.6f; particle.velocity.y *= -0.5f; }
            if (localPos.z < -0.6f) { localPos.z = -0.6f; particle.velocity.z *= -0.5f; }
            if (localPos.z > 0.6f) { localPos.z = 0.6f; particle.velocity.z *= -0.5f; }
            
            particle.position = transform.TransformPoint(localPos);
        }
    }
    
    void UpdateMeshDeformation()
    {
        // Simple vertex displacement based on nearby particle deformation
        for (int v = 0; v < currentVertices.Length; v++)
        {
            Vector3 worldVertexPos = transform.TransformPoint(originalVertices[v]);
            Vector3 totalDisplacement = Vector3.zero;
            float totalWeight = 0.0f;
            
            for (int p = 0; p < particles.Length; p++)
            {
                if (particles[p].isDeformed)
                {
                    float distance = Vector3.Distance(worldVertexPos, particles[p].position);
                    if (distance < particleRadius * 3.0f)
                    {
                        float weight = 1.0f - (distance / (particleRadius * 3.0f));
                        Vector3 displacement = (particles[p].position - particles[p].originalPosition) * 
                                             particles[p].deformationAmount * weight;
                        
                        totalDisplacement += displacement;
                        totalWeight += weight;
                    }
                }
            }
            
            if (totalWeight > 0.0f)
            {
                Vector3 worldDisplacement = totalDisplacement / totalWeight;
                Vector3 localDisplacement = transform.InverseTransformVector(worldDisplacement);
                currentVertices[v] = originalVertices[v] + localDisplacement * 0.5f;
            }
            else
            {
                // Gradually recover to original shape
                currentVertices[v] = Vector3.Lerp(currentVertices[v], originalVertices[v], 
                                                recoverySpeed * Time.deltaTime);
            }
        }
        
        // Update mesh
        clayMesh.vertices = currentVertices;
        clayMesh.RecalculateNormals();
        clayMesh.RecalculateBounds();
    }
    
    void ApplyRecovery()
    {
        for (int i = 0; i < particles.Length; i++)
        {
            if (particles[i].isDeformed)
            {
                // Gradually reduce deformation
                particles[i].deformationAmount = Mathf.Lerp(particles[i].deformationAmount, 0.0f, 
                                                          recoverySpeed * Time.deltaTime * 0.5f);
                
                // If deformation is minimal, mark as not deformed
                if (particles[i].deformationAmount < 0.01f)
                {
                    particles[i].isDeformed = false;
                    particles[i].deformationAmount = 0.0f;
                }
                
                // Slowly return to original position
                particles[i].originalPosition = Vector3.Lerp(particles[i].originalPosition, 
                    GetNearestOriginalPosition(particles[i].position), recoverySpeed * Time.deltaTime * 0.2f);
            }
        }
    }
    
    Vector3 GetNearestOriginalPosition(Vector3 currentPos)
    {
        // Find the nearest original position within the clay bounds
        Vector3 localPos = transform.InverseTransformPoint(currentPos);
        localPos = Vector3.ClampMagnitude(localPos, 0.5f);
        return transform.TransformPoint(localPos);
    }
    
    // Public methods for external control
    public void SetElasticity(float value)
    {
        elasticity = Mathf.Clamp(value, 0.1f, 2.0f);
    }
    
    public void SetPlasticity(float value)
    {
        plasticity = Mathf.Clamp01(value);
    }
    
    public void SetViscosity(float value)
    {
        viscosity = Mathf.Clamp(value, 0.1f, 5.0f);
    }
    
    public void ResetClay()
    {
        InitializeParticles();
        System.Array.Copy(originalVertices, currentVertices, originalVertices.Length);
        clayMesh.vertices = currentVertices;
        clayMesh.RecalculateNormals();
    }
    
    // Debug visualization
    void OnDrawGizmos()
    {
        if (particles != null && Application.isPlaying)
        {
            Gizmos.color = Color.blue;
            foreach (var particle in particles)
            {
                if (particle.isDeformed)
                {
                    Gizmos.color = Color.red;
                    Gizmos.DrawSphere(particle.position, particleRadius * particle.deformationAmount);
                }
                else
                {
                    Gizmos.color = Color.blue;
                    Gizmos.DrawWireSphere(particle.position, particleRadius * 0.5f);
                }
            }
        }
        
        // Draw hand interaction radius
        if (leftHand != null)
        {
            Gizmos.color = Color.green;
            Gizmos.DrawWireSphere(leftHand.position, handInteractionRadius);
        }
        
        if (rightHand != null)
        {
            Gizmos.color = Color.green;
            Gizmos.DrawWireSphere(rightHand.position, handInteractionRadius);
        }
    }
}