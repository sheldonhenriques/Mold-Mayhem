using UnityEngine;
using System.Collections.Generic;

public class SoftBody : MonoBehaviour
{
    [Header("Prefab for Visualizing Particles")]
    public GameObject particlePrefab;

    [Header("Controller Reference")]
    public Transform controllerTransform; // Drag your controller GameObject here

    [Header("Physics Settings")]
    public float gravity = -9.81f;
    public float timeStep = 0.02f;
    public int constraintIterations = 5;
    public float particleMass = 1.0f;
    public bool useGroundCollision = true;

    [Header("Interaction Modes")]
    public KeyCode spawnParticleKey = KeyCode.P; // Key to spawn individual particles
    public KeyCode connectModeKey = KeyCode.C;   // Key to toggle spring connection mode
    public KeyCode clearAllKey = KeyCode.X;      // Key to clear all particles and springs

    private List<Particle> particles = new List<Particle>();
    private List<Spring> springs = new List<Spring>();
    private bool isPinnedState = false;
    
    // New variables for spring connection mode
    private bool isConnectMode = false;
    private Particle selectedParticle1 = null;
    private Particle selectedParticle2 = null;

    void Update()
    {
        HandleKeyboardInput();
        HandleMouseClick();
        HandleControllerInput();

        foreach (var p in particles)
            p.ApplyForce(new Vector3(0, gravity, 0));

        foreach (var p in particles)
            p.UpdatePosition(timeStep);

        foreach (var p in particles)
        {
            // p.SatisfyColliderConstraint();
            p.SatisfyGroundCollision();
        }

        foreach (var s in springs)
            s.Update();
    }

    void HandleKeyboardInput()
    {
        // Spawn individual particle at random position
        if (Input.GetKeyDown(spawnParticleKey))
        {
            Vector3 randomPos = transform.position + new Vector3(
                Random.Range(-5f, 5f), 
                Random.Range(2f, 8f), 
                Random.Range(-5f, 5f)
            );
            SpawnIndividualParticle(randomPos);
        }

        // Toggle spring connection mode
        if (Input.GetKeyDown(connectModeKey))
        {
            ToggleConnectMode();
        }

        // Clear all particles and springs
        if (Input.GetKeyDown(clearAllKey))
        {
            ClearAll();
        }
    }

    void HandleMouseClick()
    {
        if (Input.GetMouseButtonDown(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out RaycastHit hit))
            {
                if (isConnectMode)
                {
                    HandleParticleSelection(hit.point);
                }
                else
                {
                    SpawnParticleAt(hit.point);
                }
            }
        }

        if (Input.GetMouseButtonDown(1))
        {
            if (isPinnedState)
            {
                foreach (var p in particles)
                    p.isPinned = false;

                if (particles.Count >= 2)
                {
                    Particle first = particles[0];
                    Particle last = particles[particles.Count - 1];
                    springs.Add(new Spring(last, first, this.transform));
                }
            }
            isPinnedState = !isPinnedState;
        }
    }

    public void SpawnParticleAt(Vector3 position)
    {
        var newParticle = new Particle(position, particleMass, particlePrefab);
        if (particles.Count == 0 || isPinnedState)
        {
            newParticle.isPinned = true;
        }
        particles.Add(newParticle);

        if (particles.Count > 1)
        {
            Particle last = particles[particles.Count - 2];
            springs.Add(new Spring(last, newParticle, this.transform));
        }
    }

    // New method to spawn particles without automatic spring connections
    public void SpawnIndividualParticle(Vector3 position)
    {
        var newParticle = new Particle(position, particleMass, particlePrefab);
        // Individual particles are not pinned by default
        newParticle.isPinned = true;
        particles.Add(newParticle);
        
        Debug.Log($"Spawned individual particle at {position}. Total particles: {particles.Count}");
    }

    void HandleParticleSelection(Vector3 clickPosition)
    {
        // Find the closest particle to the click position
        Particle closestParticle = FindClosestParticle(clickPosition);
        
        if (closestParticle == null) return;

        if (selectedParticle1 == null)
        {
            selectedParticle1 = closestParticle;
            Debug.Log("First particle selected for spring connection");
            // You could add visual feedback here (e.g., change particle color)
        }
        else if (selectedParticle2 == null && closestParticle != selectedParticle1)
        {
            selectedParticle2 = closestParticle;
            ConnectParticlesWithSpring(selectedParticle1, selectedParticle2);
            Debug.Log("Second particle selected - Spring created!");
            
            // Reset selection
            selectedParticle1 = null;
            selectedParticle2 = null;
        }
        else if (closestParticle == selectedParticle1)
        {
            // Clicked same particle, deselect
            selectedParticle1 = null;
            Debug.Log("Particle deselected");
        }
    }

    Particle FindClosestParticle(Vector3 position)
    {
        Particle closest = null;
        float minDistance = float.MaxValue;
        float maxSelectionDistance = 2f; // Maximum distance to select a particle

        foreach (var particle in particles)
        {
            float distance = Vector3.Distance(particle.position, position);
            if (distance < minDistance && distance < maxSelectionDistance)
            {
                minDistance = distance;
                closest = particle;
            }
        }

        return closest;
    }

    void ConnectParticlesWithSpring(Particle p1, Particle p2)
    {
        // Check if spring already exists between these particles
        // Note: You may need to adjust the property names based on your Spring class
        foreach (var spring in springs)
        {
            // Uncomment and adjust these lines based on your Spring class properties
            // if ((spring.particle1 == p1 && spring.particle2 == p2) ||
            //     (spring.particle1 == p2 && spring.particle2 == p1))
            // {
            //     Debug.Log("Spring already exists between these particles!");
            //     return;
            // }
        }

        springs.Add(new Spring(p1, p2, this.transform));
        Debug.Log($"Spring created between particles. Total springs: {springs.Count}");
    }

    void ToggleConnectMode()
    {
        isConnectMode = !isConnectMode;
        
        // Reset selection when toggling mode
        selectedParticle1 = null;
        selectedParticle2 = null;
        
        if (isConnectMode)
        {
            Debug.Log("Spring Connection Mode: ON - Click particles to connect them with springs");
        }
        else
        {
            Debug.Log("Normal Mode: ON - Click to spawn connected particles");
        }
    }

    void ClearAll()
    {
        // Destroy all particle GameObjects
        foreach (var particle in particles)
        {
            // Adjust this based on your Particle class structure
            // If your Particle class has a visual GameObject, uncomment and adjust:
            if (particle != null && particle.visual != null)
            {
                DestroyImmediate(particle.visual); // Destroy the visual GameObject
            }
        }

        foreach (var spring in springs)
        {
            if (spring != null && spring.line != null)
            {
                DestroyImmediate(spring.line); // Assuming spring has a visual representation
            }
        }
        
        particles.Clear();
        springs.Clear();
        selectedParticle1 = null;
        selectedParticle2 = null;
        
        Debug.Log("All particles and springs cleared");
    }

    void HandleControllerInput()
    {
        // Spawn individual particle at controller position
        if (OVRInput.GetUp(OVRInput.Button.One) && controllerTransform != null)
        {
            if (isConnectMode)
            {
                HandleParticleSelection(controllerTransform.position);
            }
            else
            {
                SpawnIndividualParticle(controllerTransform.position);
            }
        }
        
        // Toggle connection mode with controller
        if (OVRInput.GetUp(OVRInput.Button.Two))
        {
            ToggleConnectMode();
        }

        // Original connected particle spawning (moved to Button.Three)
        if (OVRInput.GetUp(OVRInput.Button.Three) && controllerTransform != null)
        {
            SpawnParticleAt(controllerTransform.position);
        }
        
        // Toggle pinned state with Button.Four
        if (OVRInput.GetUp(OVRInput.Button.Four))
        {
            if (isPinnedState)
            {
                foreach (var p in particles)
                    p.isPinned = false;

                if (particles.Count >= 2)
                {
                    Particle first = particles[0];
                    Particle last = particles[particles.Count - 1];
                    springs.Add(new Spring(last, first, this.transform));
                }
            }
            isPinnedState = !isPinnedState;
        }
    }

    void OnGUI()
    {
        GUI.Box(new Rect(10, 10, 300, 120), "SoftBody Controls");
        GUI.Label(new Rect(20, 35, 280, 20), $"P Key: Spawn Individual Particle");
        GUI.Label(new Rect(20, 55, 280, 20), $"C Key: Toggle Connect Mode ({(isConnectMode ? "ON" : "OFF")})");
        GUI.Label(new Rect(20, 75, 280, 20), $"X Key: Clear All");
        GUI.Label(new Rect(20, 95, 280, 20), $"Particles: {particles.Count}, Springs: {springs.Count}");
        
        if (isConnectMode)
        {
            GUI.Box(new Rect(10, 140, 300, 60), "Spring Connection Mode");
            GUI.Label(new Rect(20, 165, 280, 20), "Click particles to connect with springs");
            if (selectedParticle1 != null)
                GUI.Label(new Rect(20, 185, 280, 20), "First particle selected - click another");
        }
    }
}