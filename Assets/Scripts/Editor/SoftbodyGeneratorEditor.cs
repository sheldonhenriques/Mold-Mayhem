using UnityEditor;
using UnityEngine;

[CustomEditor(typeof(SoftbodyGenerator))]
public class SoftbodyGeneratorEditor : Editor
{
    public override void OnInspectorGUI()
    {
        SoftbodyGenerator softbody = (SoftbodyGenerator)target;

        softbody.debugMode = EditorGUILayout.Toggle("#Debug mod", softbody.debugMode);
        EditorGUILayout.Space();

        string[] options = new string[] { "  version 1", "  version 2" };

        softbody.gravity = EditorGUILayout.Toggle("Gravity", softbody.gravity);
        softbody.mass = EditorGUILayout.FloatField("Mass(KG)", softbody.mass);
        softbody.physicsRoughness = EditorGUILayout.FloatField("Drag (roughness)", softbody.physicsRoughness);
        softbody.softness = EditorGUILayout.FloatField("Softbody hardness", softbody.softness);
        softbody.damp = EditorGUILayout.FloatField("Softbody damper", softbody.damp);
        softbody.collissionSurfaceOffset = EditorGUILayout.FloatField("Softbody Offset", softbody.collissionSurfaceOffset);
    }
}