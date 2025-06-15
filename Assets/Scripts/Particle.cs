using UnityEngine;
public class Particle {
    public Vector3 position;
    public Vector3 previousPosition;
    public Vector3 acceleration;
    public float mass;
    public bool isPinned;
    public GameObject visual;

    public Particle(Vector3 pos, float m, GameObject prefab) {
        position = pos;
        previousPosition = pos;
        mass = m;
        acceleration = Vector3.zero;
        isPinned = false;
        visual = GameObject.Instantiate(prefab, position, Quaternion.identity);
    }

    public void ApplyForce(Vector3 force) {
        acceleration += force / mass;
    }

    public void UpdatePosition(float dt) {
        if (isPinned) return;
        Vector3 temp = position;
        position += (position - previousPosition) + acceleration * dt * dt;
        previousPosition = temp;
        acceleration = Vector3.zero;
        if (visual) visual.transform.position = position;
    }
}