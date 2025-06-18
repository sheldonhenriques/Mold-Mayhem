using UnityEngine;
public class Spring
{
    public Particle p1, p2;
    public float restLength;
    public float stiffness = 50f;
    public float damping = 0.5f;
    public LineRenderer line;

    public Spring(Particle a, Particle b, Transform parent)
    {
        p1 = a; p2 = b;
        restLength = Vector3.Distance(p1.position, p2.position);

        GameObject lineObj = new GameObject("SpringLine");
        lineObj.transform.parent = parent;
        line = lineObj.AddComponent<LineRenderer>();
        line.positionCount = 2;
        line.material = new Material(Shader.Find("Sprites/Default"));
        line.startWidth = 0.05f;
        line.endWidth = 0.05f;
        line.startColor = Color.red;
        line.endColor = Color.red;
    }

    public void ApplySpringForce()
    {
        Vector3 delta = p2.position - p1.position;
        float currentLength = delta.magnitude;
        if (currentLength == 0) return;

        Vector3 direction = delta.normalized;
        float stretch = currentLength - restLength;

        // Hooke's Law force
        Vector3 force = -stiffness * stretch * direction;

        // Approximate velocity using position differences (Verlet)
        Vector3 v1 = p1.position - p1.previousPosition;
        Vector3 v2 = p2.position - p2.previousPosition;
        Vector3 relativeVelocity = v2 - v1;

        float dampingForceMag = Vector3.Dot(relativeVelocity, direction);
        Vector3 dampingForce = -damping * dampingForceMag * direction;

        Vector3 totalForce = force + dampingForce;

        if (!p1.isPinned) p1.ApplyForce(-totalForce);
        if (!p2.isPinned) p2.ApplyForce(totalForce);
    }

    public void UpdateLine()
    {
        line.SetPosition(0, p1.position);
        line.SetPosition(1, p2.position);
    }
    
    public void Update() {
        ApplySpringForce();
        UpdateLine();
    }
}
