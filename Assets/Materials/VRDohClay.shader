Shader "Custom/VRDohClay"
{
    Properties
    {
        _BaseColor ("Base Color", Color) = (0.8, 0.4, 0.2, 1.0)
        _DeformColor ("Deform Color", Color) = (1.0, 0.6, 0.3, 1.0)
        _Roughness ("Roughness", Range(0,1)) = 0.8
        _Metallic ("Metallic", Range(0,1)) = 0.0
        _SubsurfaceColor ("Subsurface Color", Color) = (0.9, 0.5, 0.3, 1.0)
        _SubsurfaceStrength ("Subsurface Strength", Range(0,1)) = 0.3
        _NoiseScale ("Noise Scale", Range(0.1, 10)) = 2.0
        _NoiseStrength ("Noise Strength", Range(0, 0.5)) = 0.1
        _FresnelPower ("Fresnel Power", Range(0.1, 5)) = 2.0
    }
    
    SubShader
    {
        Tags { "RenderType"="Opaque" "RenderPipeline"="UniversalPipeline" }
        LOD 300
        
        Pass
        {
            Name "ForwardLit"
            Tags { "LightMode"="UniversalForward" }
            
            HLSLPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile _ _MAIN_LIGHT_SHADOWS
            #pragma multi_compile _ _MAIN_LIGHT_SHADOWS_CASCADE
            #pragma multi_compile _ _SHADOWS_SOFT
            
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Core.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Lighting.hlsl"
            
            CBUFFER_START(UnityPerMaterial)
                float4 _BaseColor;
                float4 _DeformColor;
                float4 _SubsurfaceColor;
                float _Roughness;
                float _Metallic;
                float _SubsurfaceStrength;
                float _NoiseScale;
                float _NoiseStrength;
                float _FresnelPower;
            CBUFFER_END
            
            struct Attributes
            {
                float4 positionOS : POSITION;
                float3 normalOS : NORMAL;
                float4 tangentOS : TANGENT;
                float2 uv : TEXCOORD0;
            };
            
            struct Varyings
            {
                float4 positionCS : SV_POSITION;
                float3 positionWS : TEXCOORD0;
                float3 normalWS : TEXCOORD1;
                float3 tangentWS : TEXCOORD2;
                float3 bitangentWS : TEXCOORD3;
                float2 uv : TEXCOORD4;
                float3 viewDirWS : TEXCOORD5;
                float4 shadowCoord : TEXCOORD6;
            };
            
            // Simple 3D noise function
            float noise3D(float3 pos)
            {
                return frac(sin(dot(pos, float3(12.9898, 78.233, 45.543))) * 43758.5453);
            }
            
            float fbm(float3 pos)
            {
                float value = 0.0;
                float amplitude = 0.5;
                float frequency = 1.0;
                
                for(int i = 0; i < 4; i++)
                {
                    value += amplitude * noise3D(pos * frequency);
                    amplitude *= 0.5;
                    frequency *= 2.0;
                }
                
                return value;
            }
            
            Varyings vert(Attributes input)
            {
                Varyings output;
                
                VertexPositionInputs vertexInput = GetVertexPositionInputs(input.positionOS.xyz);
                VertexNormalInputs normalInput = GetVertexNormalInputs(input.normalOS, input.tangentOS);
                
                output.positionCS = vertexInput.positionCS;
                output.positionWS = vertexInput.positionWS;
                output.normalWS = normalInput.normalWS;
                output.tangentWS = normalInput.tangentWS;
                output.bitangentWS = normalInput.bitangentWS;
                output.uv = input.uv;
                output.viewDirWS = GetWorldSpaceViewDir(vertexInput.positionWS);
                output.shadowCoord = GetShadowCoord(vertexInput);
                
                return output;
            }
            
            float4 frag(Varyings input) : SV_Target
            {
                // Base clay color with noise variation
                float3 noisePos = input.positionWS * _NoiseScale;
                float clayNoise = fbm(noisePos);
                float3 baseColor = lerp(_BaseColor.rgb, _DeformColor.rgb, clayNoise * _NoiseStrength);
                
                // Surface normal with subtle noise perturbation
                float3 normalWS = normalize(input.normalWS);
                float3 tangentWS = normalize(input.tangentWS);
                float3 bitangentWS = normalize(input.bitangentWS);
                
                // Add subtle normal variation for clay texture
                float3 noiseNormal = float3(
                    fbm(noisePos + float3(0.1, 0, 0)) - 0.5,
                    fbm(noisePos + float3(0, 0.1, 0)) - 0.5,
                    fbm(noisePos + float3(0, 0, 0.1)) - 0.5
                ) * _NoiseStrength * 0.5;
                
                normalWS = normalize(normalWS + 
                    tangentWS * noiseNormal.x + 
                    bitangentWS * noiseNormal.y + 
                    normalWS * noiseNormal.z);
                
                // View direction and lighting
                float3 viewDirWS = normalize(input.viewDirWS);
                Light mainLight = GetMainLight(input.shadowCoord);
                
                // Fresnel for rim lighting (clay-like subsurface scattering simulation)
                float fresnel = pow(1.0 - saturate(dot(normalWS, viewDirWS)), _FresnelPower);
                float3 fresnelColor = _SubsurfaceColor.rgb * fresnel * _SubsurfaceStrength;
                
                // Lambert diffuse lighting
                float NdotL = saturate(dot(normalWS, mainLight.direction));
                float3 diffuse = baseColor * mainLight.color * NdotL * mainLight.shadowAttenuation;
                
                // Simple specular reflection (very subtle for clay)
                float3 halfVector = normalize(mainLight.direction + viewDirWS);
                float NdotH = saturate(dot(normalWS, halfVector));
                float specular = pow(NdotH, lerp(1, 128, 1.0 - _Roughness)) * (1.0 - _Roughness) * 0.1;
                
                // Ambient lighting
                float3 ambient = baseColor * 0.2;
                
                // Combine all lighting
                float3 finalColor = ambient + diffuse + fresnelColor + specular;
                
                // Clay tends to be slightly desaturated
                float saturation = 0.8;
                float3 gray = dot(finalColor, float3(0.299, 0.587, 0.114));
                finalColor = lerp(gray, finalColor, saturation);
                
                return float4(finalColor, 1.0);
            }
            ENDHLSL
        }
        
        // Shadow pass
        Pass
        {
            Name "ShadowCaster"
            Tags { "LightMode"="ShadowCaster" }
            
            ZWrite On
            ZTest LEqual
            ColorMask 0
            Cull Back
            
            HLSLPROGRAM
            #pragma vertex ShadowPassVertex
            #pragma fragment ShadowPassFragment
            #include "Packages/com.unity.render-pipelines.universal/Shaders/LitInput.hlsl"
            #include "Packages/com.unity.render-pipelines.universal/Shaders/ShadowCasterPass.hlsl"
            ENDHLSL
        }
    }
    
    Fallback "Universal Render Pipeline/Lit"
}