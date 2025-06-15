using UnityEngine;
using System.Collections.Generic;
using System.Linq;

// Structure to represent a tetrahedron
[System.Serializable]
public struct Tetrahedron
{
    public int[] vertexIndices; // 4 vertex indices
    public Vector3 restCenter;
    public Matrix4x4 restTransform;
    
    public Tetrahedron(int[] indices, Vector3[] positions)
    {
        vertexIndices = indices;
        restCenter = CalculateCenter(positions, indices);
        restTransform = CalculateRestTransform(positions, indices, restCenter);
    }
    
    private static Vector3 CalculateCenter(Vector3[] positions, int[] indices)
    {
        Vector3 center = Vector3.zero;
        for (int i = 0; i < 4; i++)
        {
            center += positions[indices[i]];
        }
        return center / 4f;
    }
    
    private static Matrix4x4 CalculateRestTransform(Vector3[] positions, int[] indices, Vector3 center)
    {
        // Calculate rest-state transformation matrix for barycentric coordinates
        Vector3 v0 = positions[indices[0]] - center;
        Vector3 v1 = positions[indices[1]] - center;
        Vector3 v2 = positions[indices[2]] - center;
        Vector3 v3 = positions[indices[3]] - center;
        
        Matrix4x4 transform = new Matrix4x4();
        transform.SetRow(0, new Vector4(v0.x, v1.x, v2.x, v3.x));
        transform.SetRow(1, new Vector4(v0.y, v1.y, v2.y, v3.y));
        transform.SetRow(2, new Vector4(v0.z, v1.z, v2.z, v3.z));
        transform.SetRow(3, new Vector4(1, 1, 1, 1));
        
        return transform.inverse;
    }
}

// Structure for visual mesh vertex mapping
[System.Serializable]
public struct VertexMapping
{
    public int tetrahedronIndex;
    public Vector4 barycentricCoordinates; // 4 weights for tetrahedron vertices
}

public class TetrahedralSoftBody : MonoBehaviour
{
    [Header("Mesh References")]
    public MeshFilter visualMeshFilter;
    public MeshRenderer meshRenderer;
    
    [Header("Physics Parameters")]
    public float stiffness = 1000f;
    public float damping = 0.98f;
    public Vector3 gravity = new Vector3(0, -9.81f, 0);
    public float mass = 1f;
    
    [Header("Voxelization Settings")]
    public int voxelResolution = 8;
    public bool showTetrahedralMesh = false;
    
    // Physics simulation data
    private Vector3[] tetrahedralVertices;
    private Vector3[] tetrahedralVelocities;
    private Vector3[] tetrahedralForces;
    private Tetrahedron[] tetrahedra;
    
    // Visual mesh data
    private Mesh visualMesh;
    private Vector3[] originalVisualVertices;
    private Vector3[] deformedVisualVertices;
    private VertexMapping[] vertexMappings;
    
    // Voxel grid for tetrahedralization
    private bool[,,] voxelGrid;
    private Vector3 voxelSize;
    private Vector3 boundsMin;
    
    void Start()
    {
        InitializeMeshes();
        VoxelizeAndTetrahedralize();
        MapVisualVertices();
    }
    
    void InitializeMeshes()
    {
        visualMesh = visualMeshFilter.mesh;
        originalVisualVertices = visualMesh.vertices;
        deformedVisualVertices = new Vector3[originalVisualVertices.Length];
        System.Array.Copy(originalVisualVertices, deformedVisualVertices, originalVisualVertices.Length);
    }
    
    void VoxelizeAndTetrahedralize()
    {
        // Get mesh bounds
        Bounds bounds = visualMesh.bounds;
        boundsMin = bounds.min;
        voxelSize = bounds.size / voxelResolution;
        
        // Create voxel grid
        voxelGrid = new bool[voxelResolution, voxelResolution, voxelResolution];
        VoxelizeMesh();
        
        // Generate tetrahedral mesh from voxels
        GenerateTetrahedralMesh();
        
        // Initialize physics data
        tetrahedralVelocities = new Vector3[tetrahedralVertices.Length];
        tetrahedralForces = new Vector3[tetrahedralVertices.Length];
    }
    
    void VoxelizeMesh()
    {
        Vector3[] vertices = visualMesh.vertices;
        int[] triangles = visualMesh.triangles;
        
        // For each voxel, check if it's inside the mesh
        for (int x = 0; x < voxelResolution; x++)
        {
            for (int y = 0; y < voxelResolution; y++)
            {
                for (int z = 0; z < voxelResolution; z++)
                {
                    Vector3 voxelCenter = boundsMin + new Vector3(
                        (x + 0.5f) * voxelSize.x,
                        (y + 0.5f) * voxelSize.y,
                        (z + 0.5f) * voxelSize.z
                    );
                    
                    // Transform to local space
                    voxelCenter = transform.InverseTransformPoint(voxelCenter);
                    
                    // Use ray casting to determine if point is inside mesh
                    voxelGrid[x, y, z] = IsPointInsideMesh(voxelCenter, vertices, triangles);
                }
            }
        }
    }
    
    bool IsPointInsideMesh(Vector3 point, Vector3[] vertices, int[] triangles)
    {
        // Simple ray casting algorithm
        Vector3 rayDirection = Vector3.right;
        int intersectionCount = 0;
        
        for (int i = 0; i < triangles.Length; i += 3)
        {
            Vector3 v0 = vertices[triangles[i]];
            Vector3 v1 = vertices[triangles[i + 1]];
            Vector3 v2 = vertices[triangles[i + 2]];
            
            if (RayTriangleIntersect(point, rayDirection, v0, v1, v2))
            {
                intersectionCount++;
            }
        }
        
        return (intersectionCount % 2) == 1;
    }
    
    bool RayTriangleIntersect(Vector3 rayOrigin, Vector3 rayDirection, Vector3 v0, Vector3 v1, Vector3 v2)
    {
        // Möller-Trumbore intersection algorithm
        Vector3 edge1 = v1 - v0;
        Vector3 edge2 = v2 - v0;
        Vector3 h = Vector3.Cross(rayDirection, edge2);
        float a = Vector3.Dot(edge1, h);
        
        if (a > -0.00001f && a < 0.00001f) return false;
        
        float f = 1f / a;
        Vector3 s = rayOrigin - v0;
        float u = f * Vector3.Dot(s, h);
        
        if (u < 0f || u > 1f) return false;
        
        Vector3 q = Vector3.Cross(s, edge1);
        float v = f * Vector3.Dot(rayDirection, q);
        
        if (v < 0f || u + v > 1f) return false;
        
        float t = f * Vector3.Dot(edge2, q);
        return t > 0.00001f;
    }
    
    void GenerateTetrahedralMesh()
    {
        List<Vector3> tetVertices = new List<Vector3>();
        List<Tetrahedron> tetList = new List<Tetrahedron>();
        
        // Generate vertices from voxel grid corners
        Dictionary<Vector3Int, int> vertexMap = new Dictionary<Vector3Int, int>();
        
        for (int x = 0; x <= voxelResolution; x++)
        {
            for (int y = 0; y <= voxelResolution; y++)
            {
                for (int z = 0; z <= voxelResolution; z++)
                {
                    Vector3 pos = boundsMin + new Vector3(x * voxelSize.x, y * voxelSize.y, z * voxelSize.z);
                    Vector3Int key = new Vector3Int(x, y, z);
                    
                    vertexMap[key] = tetVertices.Count;
                    tetVertices.Add(transform.InverseTransformPoint(pos));
                }
            }
        }
        
        // Generate tetrahedra from filled voxels
        for (int x = 0; x < voxelResolution; x++)
        {
            for (int y = 0; y < voxelResolution; y++)
            {
                for (int z = 0; z < voxelResolution; z++)
                {
                    if (voxelGrid[x, y, z])
                    {
                        // Create 6 tetrahedra per cube (standard cube subdivision)
                        CreateTetrahedraFromCube(x, y, z, vertexMap, tetList, tetVertices.ToArray());
                    }
                }
            }
        }
        
        tetrahedralVertices = tetVertices.ToArray();
        tetrahedra = tetList.ToArray();
    }
    
    void CreateTetrahedraFromCube(int x, int y, int z, Dictionary<Vector3Int, int> vertexMap, 
                                  List<Tetrahedron> tetList, Vector3[] vertices)
    {
        // Get the 8 vertices of the cube
        int[] cubeVertices = new int[8];
        cubeVertices[0] = vertexMap[new Vector3Int(x, y, z)];
        cubeVertices[1] = vertexMap[new Vector3Int(x + 1, y, z)];
        cubeVertices[2] = vertexMap[new Vector3Int(x + 1, y + 1, z)];
        cubeVertices[3] = vertexMap[new Vector3Int(x, y + 1, z)];
        cubeVertices[4] = vertexMap[new Vector3Int(x, y, z + 1)];
        cubeVertices[5] = vertexMap[new Vector3Int(x + 1, y, z + 1)];
        cubeVertices[6] = vertexMap[new Vector3Int(x + 1, y + 1, z + 1)];
        cubeVertices[7] = vertexMap[new Vector3Int(x, y + 1, z + 1)];
        
        // Standard cube-to-tetrahedra subdivision (6 tetrahedra per cube)
        int[,] tetraIndices = {
            {0, 1, 2, 4}, {1, 2, 4, 5}, {2, 4, 5, 6},
            {0, 2, 3, 4}, {2, 3, 4, 7}, {2, 4, 6, 7}
        };
        
        for (int i = 0; i < 6; i++)
        {
            int[] indices = new int[4];
            for (int j = 0; j < 4; j++)
            {
                indices[j] = cubeVertices[tetraIndices[i, j]];
            }
            tetList.Add(new Tetrahedron(indices, vertices));
        }
    }
    
    void MapVisualVertices()
    {
        vertexMappings = new VertexMapping[originalVisualVertices.Length];
        
        for (int i = 0; i < originalVisualVertices.Length; i++)
        {
            Vector3 vertex = originalVisualVertices[i];
            
            // Find containing tetrahedron and compute barycentric coordinates
            for (int tetIdx = 0; tetIdx < tetrahedra.Length; tetIdx++)
            {
                Vector4 baryCoords = ComputeBarycentricCoordinates(vertex, tetIdx);
                
                // Check if point is inside tetrahedron (all barycentric coords >= 0 and sum = 1)
                if (baryCoords.x >= -0.001f && baryCoords.y >= -0.001f && 
                    baryCoords.z >= -0.001f && baryCoords.w >= -0.001f)
                {
                    vertexMappings[i] = new VertexMapping
                    {
                        tetrahedronIndex = tetIdx,
                        barycentricCoordinates = baryCoords
                    };
                    break;
                }
                
                // If not found in any tetrahedron, use closest one (fallback)
                if (tetIdx == tetrahedra.Length - 1)
                {
                    int closestTet = FindClosestTetrahedron(vertex);
                    vertexMappings[i] = new VertexMapping
                    {
                        tetrahedronIndex = closestTet,
                        barycentricCoordinates = ComputeBarycentricCoordinates(vertex, closestTet)
                    };
                }
            }
        }
    }
    
    Vector4 ComputeBarycentricCoordinates(Vector3 point, int tetrahedronIndex)
    {
        Tetrahedron tet = tetrahedra[tetrahedronIndex];
        Vector3 tetCenter = Vector3.zero;
        
        // Calculate current tetrahedron center
        for (int i = 0; i < 4; i++)
        {
            tetCenter += tetrahedralVertices[tet.vertexIndices[i]];
        }
        tetCenter /= 4f;
        
        // Transform point relative to tetrahedron center
        Vector3 localPoint = point - tetCenter;
        
        // Use the rest transform to compute barycentric coordinates
        Vector4 homogeneousPoint = new Vector4(localPoint.x, localPoint.y, localPoint.z, 1);
        Vector4 baryCoords = tet.restTransform * homogeneousPoint;
        
        return baryCoords;
    }
    
    int FindClosestTetrahedron(Vector3 point)
    {
        float minDistance = float.MaxValue;
        int closestIndex = 0;
        
        for (int i = 0; i < tetrahedra.Length; i++)
        {
            Vector3 tetCenter = Vector3.zero;
            for (int j = 0; j < 4; j++)
            {
                tetCenter += tetrahedralVertices[tetrahedra[i].vertexIndices[j]];
            }
            tetCenter /= 4f;
            
            float distance = Vector3.Distance(point, tetCenter);
            if (distance < minDistance)
            {
                minDistance = distance;
                closestIndex = i;
            }
        }
        
        return closestIndex;
    }
    
    void FixedUpdate()
    {
        SimulatePhysics();
        UpdateVisualMesh();
    }
    
    void SimulatePhysics()
    {
        float dt = Time.fixedDeltaTime;
        
        // Reset forces
        for (int i = 0; i < tetrahedralForces.Length; i++)
        {
            tetrahedralForces[i] = gravity * mass;
        }
        
        // Apply spring forces between connected vertices
        ApplyInternalForces();
        
        // Integrate
        for (int i = 0; i < tetrahedralVertices.Length; i++)
        {
            tetrahedralVelocities[i] += tetrahedralForces[i] * dt / mass;
            tetrahedralVelocities[i] *= damping;
            tetrahedralVertices[i] += tetrahedralVelocities[i] * dt;
        }
    }
    
    void ApplyInternalForces()
    {
        foreach (Tetrahedron tet in tetrahedra)
        {
            // Apply forces between all pairs of vertices in tetrahedron
            for (int i = 0; i < 4; i++)
            {
                for (int j = i + 1; j < 4; j++)
                {
                    int idx1 = tet.vertexIndices[i];
                    int idx2 = tet.vertexIndices[j];
                    
                    Vector3 pos1 = tetrahedralVertices[idx1];
                    Vector3 pos2 = tetrahedralVertices[idx2];
                    
                    Vector3 delta = pos2 - pos1;
                    float currentLength = delta.magnitude;
                    
                    // Calculate rest length (from original tetrahedron)
                    Vector3 restPos1 = transform.TransformPoint(tet.restCenter) + 
                                      transform.TransformVector(tet.restTransform.inverse.GetColumn(i));
                    Vector3 restPos2 = transform.TransformPoint(tet.restCenter) + 
                                      transform.TransformVector(tet.restTransform.inverse.GetColumn(j));
                    float restLength = Vector3.Distance(restPos1, restPos2);
                    
                    if (currentLength > 0.001f)
                    {
                        Vector3 force = delta.normalized * (currentLength - restLength) * stiffness;
                        
                        tetrahedralForces[idx1] += force;
                        tetrahedralForces[idx2] -= force;
                    }
                }
            }
        }
    }
    
    void UpdateVisualMesh()
    {
        // Deform visual mesh based on tetrahedral simulation
        for (int i = 0; i < deformedVisualVertices.Length; i++)
        {
            VertexMapping mapping = vertexMappings[i];
            Tetrahedron tet = tetrahedra[mapping.tetrahedronIndex];
            
            // Interpolate position using barycentric coordinates
            Vector3 newPosition = Vector3.zero;
            for (int j = 0; j < 4; j++)
            {
                newPosition += tetrahedralVertices[tet.vertexIndices[j]] * mapping.barycentricCoordinates[j];
            }
            
            deformedVisualVertices[i] = newPosition;
        }
        
        // Update mesh
        visualMesh.vertices = deformedVisualVertices;
        visualMesh.RecalculateNormals();
        visualMesh.RecalculateBounds();
    }
    
    public void ApplyForceAtPoint(Vector3 worldPoint, Vector3 force)
    {
        Vector3 localPoint = transform.InverseTransformPoint(worldPoint);
        
        // Find nearest tetrahedral vertex and apply force
        float minDistance = float.MaxValue;
        int nearestVertex = 0;
        
        for (int i = 0; i < tetrahedralVertices.Length; i++)
        {
            float distance = Vector3.Distance(tetrahedralVertices[i], localPoint);
            if (distance < minDistance)
            {
                minDistance = distance;
                nearestVertex = i;
            }
        }
        
        tetrahedralForces[nearestVertex] += transform.InverseTransformDirection(force);
    }
    
    void OnDrawGizmos()
    {
        if (showTetrahedralMesh && tetrahedra != null && tetrahedralVertices != null)
        {
            Gizmos.color = Color.red;
            foreach (Tetrahedron tet in tetrahedra)
            {
                // Draw tetrahedron edges
                for (int i = 0; i < 4; i++)
                {
                    for (int j = i + 1; j < 4; j++)
                    {
                        Vector3 pos1 = transform.TransformPoint(tetrahedralVertices[tet.vertexIndices[i]]);
                        Vector3 pos2 = transform.TransformPoint(tetrahedralVertices[tet.vertexIndices[j]]);
                        Gizmos.DrawLine(pos1, pos2);
                    }
                }
            }
        }
    }
}