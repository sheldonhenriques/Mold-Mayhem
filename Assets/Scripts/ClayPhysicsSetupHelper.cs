using UnityEngine;
using UnityEditor;

public class ClayPhysicsSetupHelper : MonoBehaviour
{
    [Header("Setup Options")]
    public bool useOptimizedVersion = true;
    public bool autoFindVRHands = true;
    public Material clayMaterial;
    
    [Header("Clay Presets")]
    public ClayPreset clayPreset = ClayPreset.Soft;
    
    public enum ClayPreset
    {
        VerySoft,
        Soft,
        Medium,
        Firm,
        Hard
    }
    
    [System.Serializable]
    public struct ClayProperties
    {
        public float elasticity;
        public float plasticity;
        public float viscosity;
        public float particleRadius;
        public int particleCount;
        public float maxDeformation;
        public float recoverySpeed;
    }
    
    private ClayProperties[] presets = new ClayProperties[]
    {
        // Very Soft
        new ClayProperties { elasticity = 0.3f, plasticity = 0.8f, viscosity = 0.5f, particleRadius = 0.08f, particleCount = 300, maxDeformation = 0.8f, recoverySpeed = 1.0f },
        // Soft
        new ClayProperties { elasticity = 0.5f, plasticity = 0.6f, viscosity = 1.0f, particleRadius = 0.06f, particleCount = 400, maxDeformation = 0.6f, recoverySpeed = 1.5f },
        // Medium
        new ClayProperties { elasticity = 0.8f, plasticity = 0.4f, viscosity = 1.5f, particleRadius = 0.05f, particleCount = 500, maxDeformation = 0.4f, recoverySpeed = 2.0f },
        // Firm
        new ClayProperties { elasticity = 1.2f, plasticity = 0.2f, viscosity = 2.0f, particleRadius = 0.04f, particleCount = 600, maxDeformation = 0.3f, recoverySpeed = 3.0f },
        // Hard
        new ClayProperties { elasticity = 1.8f, plasticity = 0.1f, viscosity = 3.0f, particleRadius = 0.03f, particleCount = 800, maxDeformation = 0.2f, recoverySpeed = 4.0f }
    };
    
    [ContextMenu("Setup Clay Physics")]
    public void SetupClayPhysics()
    {
        // Remove existing clay physics components
        var existingClay = GetComponent<ClayPhysics>();
        var existingOptimizedClay = GetComponent<OptimizedClayPhysics>();
        
        if (existingClay != null)
        {
            DestroyImmediate(existingClay);
        }
        
        if (existingOptimizedClay != null)
        {
            DestroyImmediate(existingOptimizedClay);
        }
        
        // Ensure we have required components
        if (GetComponent<MeshFilter>() == null)
        {
            gameObject.AddComponent<MeshFilter>();
        }
        
        if (GetComponent<MeshRenderer>() == null)
        {
            gameObject.AddComponent<MeshRenderer>();
        }
        
        // Create a cube mesh if none exists
        MeshFilter meshFilter = GetComponent<MeshFilter>();
        if (meshFilter.sharedMesh == null)
        {
            GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
            meshFilter.sharedMesh = cube.GetComponent<MeshFilter>().sharedMesh;
            DestroyImmediate(cube);
        }
        
        // Apply clay material
        MeshRenderer renderer = GetComponent<MeshRenderer>();
        if (clayMaterial != null)
        {
            renderer.material = clayMaterial;
        }
        else
        {
            // Create a basic clay material
            Material defaultClay = new Material(Shader.Find("Standard"));
            defaultClay.name = "Clay Material";
            defaultClay.color = new Color(0.8f, 0.6f, 0.4f); // Clay color
            defaultClay.SetFloat("_Metallic", 0.0f);
            defaultClay.SetFloat("_Smoothness", 0.3f);
            renderer.material = defaultClay;
        }
        
        // Add clay physics component
        if (useOptimizedVersion)
        {
            OptimizedClayPhysics clayPhysics = gameObject.AddComponent<OptimizedClayPhysics>();
            ApplyPreset(clayPhysics);
            
            // Try to find and assign compute shader
            ComputeShader computeShader = Resources.Load<ComputeShader>("ClayPhysics");
            if (computeShader != null)
            {
                clayPhysics.clayComputeShader = computeShader;
            }
            else
            {
                Debug.LogWarning("ClayPhysics compute shader not found in Resources folder. Using CPU fallback.");
            }
        }
        else
        {
            ClayPhysics clayPhysics = gameObject.AddComponent<ClayPhysics>();
            ApplyBasicPreset(clayPhysics);
        }
        
        // Auto-find VR hands
        if (autoFindVRHands)
        {
            FindAndAssignVRHands();
        }
        
        Debug.Log($"Clay physics setup complete with {clayPreset} preset!");
    }
    
    void ApplyPreset(OptimizedClayPhysics clayPhysics)
    {
        ClayProperties preset = presets[(int)clayPreset];
        
        clayPhysics.elasticity = preset.elasticity;
        clayPhysics.plasticity = preset.plasticity;
        clayPhysics.viscosity = preset.viscosity;
        clayPhysics.particleRadius = preset.particleRadius;
        clayPhysics.particleCount = preset.particleCount;
        clayPhysics.maxDeformation = preset.maxDeformation;
        clayPhysics.recoverySpeed = preset.recoverySpeed;
        clayPhysics.handInteractionRadius = 0.1f;
    }
    
    void ApplyBasicPreset(ClayPhysics clayPhysics)
    {
        ClayProperties preset = presets[(int)clayPreset];
        
        clayPhysics.elasticity = preset.elasticity;
        clayPhysics.plasticity = preset.plasticity;
        clayPhysics.viscosity = preset.viscosity;
        clayPhysics.particleRadius = preset.particleRadius;
        clayPhysics.particleCount = Mathf.Min(preset.particleCount, 300); // Limit for CPU version
        clayPhysics.maxDeformation = preset.maxDeformation;
        clayPhysics.recoverySpeed = preset.recoverySpeed;
        clayPhysics.handInteractionRadius = 0.1f;
    }
    
    void FindAndAssignVRHands()
    {
        // Common VR hand object names
        string[] leftHandNames = { "LeftHand", "Left Hand", "HandLeft", "OVRHandPrefab_Left", "LeftControllerAnchor" };
        string[] rightHandNames = { "RightHand", "Right Hand", "HandRight", "OVRHandPrefab_Right", "RightControllerAnchor" };
        
        Transform leftHand = null;
        Transform rightHand = null;
        
        // Search for left hand
        foreach (string name in leftHandNames)
        {
            GameObject leftHandObj = GameObject.Find(name);
            if (leftHandObj != null)
            {
                leftHand = leftHandObj.transform;
                break;
            }
        }
        
        // Search for right hand
        foreach (string name in rightHandNames)
        {
            GameObject rightHandObj = GameObject.Find(name);
            if (rightHandObj != null)
            {
                rightHand = rightHandObj.transform;
                break;
            }
        }
        
        // Assign to clay physics
        if (useOptimizedVersion)
        {
            OptimizedClayPhysics clayPhysics = GetComponent<OptimizedClayPhysics>();
            if (clayPhysics != null)
            {
                clayPhysics.leftHand = leftHand;
                clayPhysics.rightHand = rightHand;
            }
        }
        else
        {
            ClayPhysics clayPhysics = GetComponent<ClayPhysics>();
            if (clayPhysics != null)
            {
                clayPhysics.leftHand = leftHand;
                clayPhysics.rightHand = rightHand;
            }
        }
        
        if (leftHand != null || rightHand != null)
        {
            Debug.Log($"VR Hands found and assigned: Left={leftHand?.name}, Right={rightHand?.name}");
        }
        else
        {
            Debug.LogWarning("No VR hands found. Please assign them manually in the inspector.");
        }
    }
    
    [ContextMenu("Create Clay Material")]
    public void CreateClayMaterial()
    {
        Material clayMat = new Material(Shader.Find("Standard"));
        clayMat.name = "Clay Material";
        
        // Clay-like properties
        clayMat.color = new Color(0.8f, 0.6f, 0.4f);
        clayMat.SetFloat("_Metallic", 0.0f);
        clayMat.SetFloat("_Smoothness", 0.3f);
        clayMat.SetFloat("_BumpScale", 1.0f);
        
        // Save as asset
        #if UNITY_EDITOR
        string path = "Assets/ClayMaterial.mat";
        AssetDatabase.CreateAsset(clayMat, path);
        AssetDatabase.SaveAssets();
        AssetDatabase.Refresh();
        
        clayMaterial = clayMat;
        Debug.Log($"Clay material created at {path}");
        #endif
    }
    
    [ContextMenu("Test Clay Properties")]
    public void TestClayProperties()
    {
        if (useOptimizedVersion)
        {
            OptimizedClayPhysics clayPhysics = GetComponent<OptimizedClayPhysics>();
            if (clayPhysics != null)
            {
                Debug.Log($"Clay Properties - Elasticity: {clayPhysics.elasticity}, Plasticity: {clayPhysics.plasticity}, Particles: {clayPhysics.particleCount}");
                Debug.Log($"Active Particles: {clayPhysics.GetActiveParticleCount()}, Performance: {clayPhysics.GetPerformanceMetric():F1} FPS");
            }
        }
        else
        {
            ClayPhysics clayPhysics = GetComponent<ClayPhysics>();
            if (clayPhysics != null)
            {
                Debug.Log($"Clay Properties - Elasticity: {clayPhysics.elasticity}, Plasticity: {clayPhysics.plasticity}");
            }
        }
    }
}

#if UNITY_EDITOR
[CustomEditor(typeof(ClayPhysicsSetupHelper))]
public class ClayPhysicsSetupHelperEditor : Editor
{
    public override void OnInspectorGUI()
    {
        DrawDefaultInspector();
        
        EditorGUILayout.Space();
        
        ClayPhysicsSetupHelper helper = (ClayPhysicsSetupHelper)target;
        
        if (GUILayout.Button("Setup Clay Physics", GUILayout.Height(30)))
        {
            helper.SetupClayPhysics();
        }
        
        EditorGUILayout.Space();
        
        EditorGUILayout.BeginHorizontal();
        if (GUILayout.Button("Create Clay Material"))
        {
            helper.CreateClayMaterial();
        }
        
        if (GUILayout.Button("Test Properties"))
        {
            helper.TestClayProperties();
        }
        EditorGUILayout.EndHorizontal();
        
        EditorGUILayout.Space();
        
        EditorGUILayout.HelpBox(
            "1. Attach this script to a cube GameObject\n" +
            "2. Choose your clay preset and options\n" +
            "3. Click 'Setup Clay Physics'\n" +
            "4. Assign VR hand references if not auto-detected\n" +
            "5. For optimized version, place ClayPhysics.compute in Resources folder",
            MessageType.Info
        );
    }
}
#endif