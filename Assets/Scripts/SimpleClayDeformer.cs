using UnityEngine;

public class SimpleClayDeformer : MonoBehaviour
{
    [Header("Deformation Settings")]
    public float deformationStrength = 0.1f;
    public float deformationRadius = 1f;
    
    private Mesh originalMesh;
    private Mesh workingMesh;
    private Vector3[] originalVertices;
    private Vector3[] modifiedVertices;
    
    void Start()
    {
        MeshFilter meshFilter = GetComponent<MeshFilter>();
        originalMesh = meshFilter.mesh;
        
        // Create a copy to work with
        workingMesh = Instantiate(originalMesh);
        meshFilter.mesh = workingMesh;
        
        originalVertices = originalMesh.vertices;
        modifiedVertices = new Vector3[originalVertices.Length];
        System.Array.Copy(originalVertices, modifiedVertices, originalVertices.Length);
    }
    
    void Update()
    {
        // Simple mouse interaction for testing
        if (Input.GetMouseButton(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            RaycastHit hit;
            
            if (Physics.Raycast(ray, out hit) && hit.collider.gameObject == gameObject)
            {
                DeformMesh(hit.point);
            }
        }
        
        if (Input.GetKeyDown(KeyCode.R))
        {
            ResetMesh();
        }
    }
    
    void DeformMesh(Vector3 hitPoint)
    {
        Vector3 localHitPoint = transform.InverseTransformPoint(hitPoint);
        
        for (int i = 0; i < modifiedVertices.Length; i++)
        {
            float distance = Vector3.Distance(modifiedVertices[i], localHitPoint);
            
            if (distance < deformationRadius)
            {
                float strength = 1f - (distance / deformationRadius);
                Vector3 deformation = (localHitPoint - modifiedVertices[i]).normalized * deformationStrength * strength;
                modifiedVertices[i] += deformation;
            }
        }
        
        workingMesh.vertices = modifiedVertices;
        workingMesh.RecalculateNormals();
        workingMesh.RecalculateBounds();
    }
    
    void ResetMesh()
    {
        System.Array.Copy(originalVertices, modifiedVertices, originalVertices.Length);
        workingMesh.vertices = modifiedVertices;
        workingMesh.RecalculateNormals();
        workingMesh.RecalculateBounds();
    }
}