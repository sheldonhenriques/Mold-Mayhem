using UnityEngine;
public class Particle
{
    public Vector3 position;
    public Vector3 previousPosition;
    public Vector3 acceleration;
    public float mass;
    public bool isPinned;
    public GameObject visual;
    public float radius = 0.2f;

    public float damping = 0.98f; // Global damping factor
    public float restitution = 0.1f; // Bounce factor (0 = no bounce, 1 = perfect bounce)

    public Particle(Vector3 pos, float m, GameObject prefab)
    {
        position = pos;
        previousPosition = pos;
        mass = m;
        acceleration = Vector3.zero;
        isPinned = false;
        visual = GameObject.Instantiate(prefab, position, Quaternion.identity);
    }

    public void ApplyForce(Vector3 force)
    {
        acceleration += force / mass;
    }

    public void UpdatePosition(float dt)
    {
        if (isPinned) return;
        Vector3 velocity = (position - previousPosition) * damping;
        Vector3 temp = position;
        position += velocity + acceleration * dt * dt;
        previousPosition = temp;
        acceleration = Vector3.zero;
        if (visual) visual.transform.position = position;
    }
    
    public void SatisfyGroundCollision(float groundHeight = 0f, float restitution = 0.2f, float friction = 0.8f)
    {
        if (isPinned) return;
        
        if (position.y - radius < groundHeight)
        {
            position.y = groundHeight + radius;
            
            // Get velocity from Verlet integration
            Vector3 velocity = position - previousPosition;
            
            // Only apply bounce if moving downward and with significant speed
            if (velocity.y < -0.005f) // Small threshold to prevent micro-bounces
            {
                // Separate velocity into normal (vertical) and tangential (horizontal) components
                Vector3 normalVelocity = new Vector3(0, velocity.y, 0);
                Vector3 tangentialVelocity = new Vector3(velocity.x, 0, velocity.z);
                
                // Apply restitution to normal velocity (bounce)
                normalVelocity *= -restitution;
                
                // Apply friction to tangential velocity (sliding)
                tangentialVelocity *= (1f - friction);
                
                // Combine and update previous position
                Vector3 newVelocity = normalVelocity + tangentialVelocity;
                previousPosition = position - newVelocity;
            }
            else
            {
                // If barely moving vertically, just stop the vertical motion
                velocity.y = 0;
                previousPosition = position - velocity;
            }
        }
    }
}