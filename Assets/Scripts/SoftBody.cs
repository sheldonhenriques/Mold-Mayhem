using UnityEngine;
using System.Collections.Generic;

public class SoftBody : MonoBehaviour
{
    [Header("Prefab for Visualizing Particles")]
    public GameObject particlePrefab;

    [Header("Controller Reference")]
    public Transform controllerTransform; // Drag your controller GameObject here

    public float gravity = -9.81f;
    public float timeStep = 0.02f;
    public int constraintIterations = 5;
    public float particleMass = 1.0f;

    private List<Particle> particles = new List<Particle>();
    private List<Spring> springs = new List<Spring>();

    void Update()
    {
        HandleMouseClick();
        HandleControllerInput();

        foreach (var p in particles)
            p.ApplyForce(new Vector3(0, gravity, 0));

        foreach (var p in particles)
            p.UpdatePosition(timeStep);

        for (int i = 0; i < constraintIterations; i++)
            foreach (var s in springs)
                s.SolveConstraint();

        foreach (var s in springs)
            s.Update();
    }

    void HandleMouseClick()
    {
        if (Input.GetMouseButtonDown(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out RaycastHit hit))
            {
                SpawnParticleAt(hit.point);
            }
        }
    }

    public void SpawnParticleAt(Vector3 position)
    {
        var newParticle = new Particle(position, particleMass, particlePrefab);
        if (particles.Count == 0)
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
    

    void HandleControllerInput()
    {
        if (OVRInput.GetUp(OVRInput.Button.One) && controllerTransform != null)
        {
            SpawnParticleAt(controllerTransform.position);
        }
    }
}
