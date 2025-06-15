using UnityEngine;
using System.Collections.Generic;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class OptimizedClayPhysics : MonoBehaviour
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
    
    [Range(100, 2000)]
    public int particleCount = 500;
    
    [Header("Deformation Settings")]
    public float maxDeformation = 0.5f;
    public float recoverySpeed = 2.0f;
    
    [Header("VR Hand References")]
    public Transform leftHand;
    public Transform rightHand;
    public float handInteractionRadius = 0.1f;
    
    [Header("Compute Shader")]
    public ComputeShader clayComputeShader;
    
    // Compute shader kernels
    private int updateParticlesKernel;
    private int applyConstraintsKernel;
    
    // Compute buffers
    private ComputeBuffer particleBuffer;
    private ClayParticleData[] particleData;
    
    // Mesh data
    private Vector3[] originalVertices;
    private Vector3[] currentVertices;
    private MeshFilter meshFilter;
    private Mesh clayMesh;
    
    // VR hand tracking
    private Vector3 leftHandVelocity;
    private Vector3 rightHandVelocity;
    private Vector3 lastLeftHandPos;
    private Vector3 lastRightHandPos;
    
    // Performance optimization
    private float lastUpdateTime;
    private const float UPDATE_INTERVAL = 0.016f; // ~60fps
    
    [System.Serializable]
    [System.Runtime.InteropServices.StructLayout(System.Runtime.InteropServices.LayoutKind.Sequential)]
    public struct ClayParticleData
    {
        public Vector3 position;
        public Vector3 velocity;
        public Vector3 force;
        public float mass;
        public float pressure;
        public float density;
        public int isDeformed;
        public float deformationAmount;
        public Vector3 originalPosition;
    }
    
    void Start()
    {
        InitializeClaySystem();
        InitializeComputeShader();
    }
    
    void InitializeClaySystem()
    {
        meshFilter = GetComponent<MeshFilter>();
        clayMesh = Instantiate(meshFilter.sharedMesh);
        meshFilter.mesh = clayMesh;
        
        originalVertices = clayMesh.vertices;
        currentVertices = new Vector3[originalVertices.Length];
        System.Array.Copy(originalVertices, currentVertices, originalVertices.Length);
        
        InitializeParticles();
        
        // Initialize hand tracking
        if (leftHand != null) lastLeftHandPos = leftHand.position;
        if (rightHand != null) lastRightHandPos = rightHand.position;
        
        Debug.Log($"Optimized Clay Physics initialized with {particleCount} particles");
    }
    
    void InitializeParticles()
    {
        particleData = new ClayParticleData[particleCount];
        
        for (int i = 0; i < particleCount; i++)
        {
            Vector3 randomPos = GetRandomPositionInBounds();
            particleData[i] = new ClayParticleData
            {
                position = randomPos,
                velocity = Vector3.zero,
                force = Vector3.zero,
                mass = 1.0f,
                pressure = 0.0f,
                density = 1.0f,
                isDeformed = 0,
                deformationAmount = 0.0f,
                originalPosition = randomPos
            };
        }
    }
    
    Vector3 GetRandomPositionInBounds()
    {
        Vector3 localPos = new Vector3(
            Random.Range(-0.4f, 0.4f),
            Random.Range(-0.4f, 0.4f),
            Random.Range(-0.4f, 0.4f)
        );
        
        return transform.TransformPoint(localPos);
    }
    
    void InitializeComputeShader()
    {
        if (clayComputeShader == null)
        {
            Debug.LogError("Clay Compute Shader not assigned! Falling back to CPU simulation.");
            return;
        }
        
        // Find kernels
        updateParticlesKernel = clayComputeShader.FindKernel("UpdateParticles");
        applyConstraintsKernel = clayComputeShader.FindKernel("ApplyConstraints");
        
        // Create compute buffer
        particleBuffer = new ComputeBuffer(particleCount, System.Runtime.InteropServices.Marshal.SizeOf<ClayParticleData>());
        particleBuffer.SetData(particleData);
        
        // Set buffer to compute shader
        clayComputeShader.SetBuffer(updateParticlesKernel, "particles", particleBuffer);
        clayComputeShader.SetBuffer(applyConstraintsKernel, "particles", particleBuffer);
        
        // Set constant properties
        clayComputeShader.SetFloat("particleRadius", particleRadius);
        clayComputeShader.SetFloat("maxDeformation", maxDeformation);
        clayComputeShader.SetInt("particleCount", particleCount);
        clayComputeShader.SetFloat("handInteractionRadius", handInteractionRadius);
        
        // Set transform matrices
        clayComputeShader.SetMatrix("worldToLocal", transform.worldToLocalMatrix);
        clayComputeShader.SetMatrix("localToWorld", transform.localToWorldMatrix);
    }
    
    void Update()
    {
        if (Time.time - lastUpdateTime >= UPDATE_INTERVAL)
        {
            UpdateHandTracking();
            
            if (clayComputeShader != null && particleBuffer != null)
            {
                UpdateClayPhysicsGPU();
            }
            else
            {
                UpdateClayPhysicsCPU();
            }
            
            UpdateMeshDeformation();
            lastUpdateTime = Time.time;
        }
    }
    
    void UpdateHandTracking()
    {
        // Calculate hand velocities
        if (leftHand != null)
        {
            leftHandVelocity = (leftHand.position - lastLeftHandPos) / Time.deltaTime;
            lastLeftHandPos = leftHand.position;
        }
        
        if (rightHand != null)
        {
            rightHandVelocity = (rightHand.position - lastRightHandPos) / Time.deltaTime;
            lastRightHandPos = rightHand.position;
        }
    }
    
    void UpdateClayPhysicsGPU()
    {
        // Update shader properties
        clayComputeShader.SetFloat("elasticity", elasticity);
        clayComputeShader.SetFloat("plasticity", plasticity);
        clayComputeShader.SetFloat("viscosity", viscosity);
        clayComputeShader.SetFloat("deltaTime", Time.deltaTime);
        
        // Update hand positions
        if (leftHand != null)
        {
            clayComputeShader.SetVector("leftHandPos", leftHand.position);
            clayComputeShader.SetVector("leftHandVelocity", leftHandVelocity);
        }
        else
        {
            clayComputeShader.SetVector("leftHandPos", Vector3.zero);
            clayComputeShader.SetVector("leftHandVelocity", Vector3.zero);
        }
        
        if (rightHand != null)
        {
            clayComputeShader.SetVector("rightHandPos", rightHand.position);
            clayComputeShader.SetVector("rightHandVelocity", rightHandVelocity);
        }
        else
        {
            clayComputeShader.SetVector("rightHandPos", Vector3.zero);
            clayComputeShader.SetVector("rightHandVelocity", Vector3.zero);
        }
        
        // Update transform matrices
        clayComputeShader.SetMatrix("worldToLocal", transform.worldToLocalMatrix);
        clayComputeShader.SetMatrix("localToWorld", transform.localToWorldMatrix);
        
        // Dispatch compute shaders
        int threadGroups = Mathf.CeilToInt(particleCount / 64.0f);
        
        // Apply constraints first
        clayComputeShader.Dispatch(applyConstraintsKernel, threadGroups, 1, 1);
        
        // Then update particle positions
        clayComputeShader.Dispatch(updateParticlesKernel, threadGroups, 1, 1);
        
        // Read back data for mesh deformation
        particleBuffer.GetData(particleData);
    }
    
    void UpdateClayPhysicsCPU()
    {
        // Fallback CPU implementation
        // Reset forces
        for (int i = 0; i < particleData.Length; i++)
        {
            var particle = particleData[i];
            particle.force = Vector3.zero;
            particleData[i] = particle;
        }
        
        // Apply elastoplasticity constraints
        ApplyElastoplasticityConstraintsCPU();
        
        // Handle VR hand interactions
        HandleVRHandInteractionsCPU();
        
        // Update particle positions
        UpdateParticlePositionsCPU();
        
        // Apply recovery
        ApplyRecoveryCPU();
    }
    
    void ApplyElastoplasticityConstraintsCPU()
    {
        for (int i = 0; i < particleData.Length; i++)
        {
            for (int j = i + 1; j < particleData.Length; j++)
            {
                var particleA = particleData[i];
                var particleB = particleData[j];
                
                Vector3 displacement = particleB.position - particleA.position;
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
                        particleA.isDeformed = 1;
                        particleB.isDeformed = 1;
                        particleA.deformationAmount = Mathf.Min(strain, maxDeformation);
                        particleB.deformationAmount = Mathf.Min(strain, maxDeformation);
                        
                        // Reduce elastic force for plastic behavior
                        elasticForce *= (1.0f - plasticity);
                    }
                    
                    // Apply viscosity damping
                    Vector3 relativeVelocity = particleB.velocity - particleA.velocity;
                    Vector3 viscousForce = relativeVelocity * viscosity * -1.0f;
                    
                    // Apply forces
                    particleA.force -= elasticForce + viscousForce;
                    particleB.force += elasticForce + viscousForce;
                    
                    particleData[i] = particleA;
                    particleData[j] = particleB;
                }
            }
        }
    }
    
    void HandleVRHandInteractionsCPU()
    {
        HandleHandInteractionCPU(leftHand, leftHandVelocity);
        HandleHandInteractionCPU(rightHand, rightHandVelocity);
    }
    
    void HandleHandInteractionCPU(Transform hand, Vector3 handVelocity)
    {
        if (hand == null) return;
        
        Vector3 handWorldPos = hand.position;
        
        for (int i = 0; i < particleData.Length; i++)
        {
            var particle = particleData[i];
            float distance = Vector3.Distance(particle.position, handWorldPos);
            
            if (distance < handInteractionRadius)
            {
                // Calculate interaction force
                Vector3 direction = (particle.position - handWorldPos).normalized;
                float forceMagnitude = (handInteractionRadius - distance) / handInteractionRadius;
                
                // Apply deformation
                particle.isDeformed = 1;
                particle.deformationAmount = Mathf.Max(particle.deformationAmount, forceMagnitude * maxDeformation);
                
                // Push particles away from hand
                Vector3 pushForce = direction * forceMagnitude * 10.0f;
                particle.force += pushForce;
                
                // Add hand velocity influence
                particle.velocity += handVelocity * 0.1f;
                
                particleData[i] = particle;
            }
        }
    }
    
    void UpdateParticlePositionsCPU()
    {
        float deltaTime = Time.deltaTime;
        
        for (int i = 0; i < particleData.Length; i++)
        {
            var particle = particleData[i];
            
            // Integrate velocity
            particle.velocity += (particle.force / particle.mass) * deltaTime;
            
            // Apply air resistance
            particle.velocity *= 0.98f;
            
            // Integrate position
            particle.position += particle.velocity * deltaTime;
            
            // Keep particles within bounds
            ConstrainToBoundsCPU(ref particle);
            
            particleData[i] = particle;
        }
    }
    
    void ConstrainToBoundsCPU(ref ClayParticleData particle)
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
    
    void ApplyRecoveryCPU()
    {
        for (int i = 0; i < particleData.Length; i++)
        {
            var particle = particleData[i];
            
            if (particle.isDeformed == 1)
            {
                // Gradually reduce deformation
                particle.deformationAmount = Mathf.Lerp(particle.deformationAmount, 0.0f, 
                                                      recoverySpeed * Time.deltaTime * 0.5f);
                
                // If deformation is minimal, mark as not deformed
                if (particle.deformationAmount < 0.01f)
                {
                    particle.isDeformed = 0;
                    particle.deformationAmount = 0.0f;
                }
                
                particleData[i] = particle;
            }
        }
    }
    
    void UpdateMeshDeformation()
    {
        // Advanced mesh deformation using particle influence
        for (int v = 0; v < currentVertices.Length; v++)
        {
            Vector3 worldVertexPos = transform.TransformPoint(originalVertices[v]);
            Vector3 totalDisplacement = Vector3.zero;
            float totalWeight = 0.0f;
            
            // Find nearby deformed particles
            for (int p = 0; p < particleData.Length; p++)
            {
                if (particleData[p].isDeformed == 1)
                {
                    float distance = Vector3.Distance(worldVertexPos, particleData[p].position);
                    float influenceRadius = particleRadius * 4.0f;
                    
                    if (distance < influenceRadius)
                    {
                        float weight = 1.0f - (distance / influenceRadius);
                        weight = weight * weight; // Smooth falloff
                        
                        Vector3 displacement = (particleData[p].position - particleData[p].originalPosition) * 
                                             particleData[p].deformationAmount * weight;
                        
                        totalDisplacement += displacement;
                        totalWeight += weight;
                    }
                }
            }
            
            if (totalWeight > 0.0f)
            {
                Vector3 worldDisplacement = totalDisplacement / totalWeight;
                Vector3 localDisplacement = transform.InverseTransformVector(worldDisplacement);
                
                // Apply displacement with smoothing
                Vector3 targetVertex = originalVertices[v] + localDisplacement * 0.3f;
                currentVertices[v] = Vector3.Lerp(currentVertices[v], targetVertex, Time.deltaTime * 5.0f);
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
    
    // Public methods for runtime control
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
    
    public void SetHandInteractionRadius(float radius)
    {
        handInteractionRadius = radius;
        if (clayComputeShader != null)
        {
            clayComputeShader.SetFloat("handInteractionRadius", handInteractionRadius);
        }
    }
    
    public void ResetClay()
    {
        InitializeParticles();
        
        if (particleBuffer != null)
        {
            particleBuffer.SetData(particleData);
        }
        
        System.Array.Copy(originalVertices, currentVertices, originalVertices.Length);
        clayMesh.vertices = currentVertices;
        clayMesh.RecalculateNormals();
        
        Debug.Log("Clay physics reset");
    }
    
    // VR-specific methods
    public void SetVRHands(Transform leftHandTransform, Transform rightHandTransform)
    {
        leftHand = leftHandTransform;
        rightHand = rightHandTransform;
        
        if (leftHand != null) lastLeftHandPos = leftHand.position;
        if (rightHand != null) lastRightHandPos = rightHand.position;
        
        Debug.Log("VR hands assigned to clay physics");
    }
    
    public void EnableHapticFeedback(float intensity = 0.1f)
    {
        // This would integrate with Meta SDK for haptic feedback
        // when hands interact with clay
        StartCoroutine(HapticFeedbackCoroutine(intensity));
    }
    
    private System.Collections.IEnumerator HapticFeedbackCoroutine(float intensity)
    {
        // Placeholder for Meta SDK haptic integration
        // You would implement actual haptic feedback here
        yield return new WaitForSeconds(0.1f);
    }
    
    // Performance monitoring
    public int GetActiveParticleCount()
    {
        int activeCount = 0;
        foreach (var particle in particleData)
        {
            if (particle.isDeformed == 1) activeCount++;
        }
        return activeCount;
    }
    
    public float GetPerformanceMetric()
    {
        return 1.0f / Time.deltaTime; // Simple FPS metric
    }
    
    // Cleanup
    void OnDestroy()
    {
        if (particleBuffer != null)
        {
            particleBuffer.Release();
            particleBuffer = null;
        }
    }
    
    // Debug visualization
    void OnDrawGizmos()
    {
        if (particleData != null && Application.isPlaying)
        {
            int deformedCount = 0;
            
            foreach (var particle in particleData)
            {
                if (particle.isDeformed == 1)
                {
                    Gizmos.color = Color.Lerp(Color.yellow, Color.red, particle.deformationAmount);
                    Gizmos.DrawSphere(particle.position, particleRadius * (1.0f + particle.deformationAmount));
                    deformedCount++;
                }
                else
                {
                    Gizmos.color = Color.blue * 0.3f;
                    Gizmos.DrawWireSphere(particle.position, particleRadius * 0.5f);
                }
            }
        }
        
        // Draw hand interaction zones
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