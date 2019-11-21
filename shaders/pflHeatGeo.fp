#version 130 
in vec3 viewDir;

in vec3 camSpacePosition;
in vec3 camSpaceNormal;
in float utv;
in float phiv;

uniform vec3 lightCamSpacePosition;
uniform vec3 lightColor;
uniform vec3 materialDiffuse;
uniform int shininess;
uniform float power;
//uniform float maxphi;

void main(void) {
	
	vec3 v = normalize(-camSpacePosition);
    // normal of fragment in camera space
	vec3 n = normalize(camSpaceNormal);
	vec3 l = normalize(lightCamSpacePosition - camSpacePosition);

	// TODO: 4, 11 Implement your GLSL per fragement lighting, heat colouring, and distance stripes here!


    vec3 reflectDir = reflect(l , n.xyz);
    float spec = max(dot(viewDir, reflectDir), 0.0);

    float diff = max(dot(n, l), 0.0);

    if (diff == 0.0)
        spec = 0.0;
    else
        spec = pow(spec, shininess);

    float spec_strength = 0.1;
    vec3 specular = spec_strength * spec * lightColor;

    vec3 diffuse = diff * materialDiffuse  * lightColor;

    float amb_strength = 0.1;
    vec3 ambient = amb_strength * lightColor;
    vec3 light = (diffuse + specular + ambient);










    float heat = pow(utv, power);





    if (phiv != 0){

        // the middle of each line is at phiv % middle_line_dist
        float middle_line_dist = 4.0;
        float full_range = 2.0;
        float phi_mod = mod(phiv*40.0, middle_line_dist + full_range/2.0);

        if (phi_mod > middle_line_dist - full_range/2.0){

            float add = phi_mod - 3.0;
            // range I add values on is: middle-1 -> middle+1
            // (div by range = 2)
            // add linearly increases from 0 -> 1

            add = add/2.0;
            float step = smoothstep(0.22, 0.30, add) * (1 - smoothstep(0.62, 0.70, add));



            // having it subtract by heat/2 means close points will be more red
            // gl_FragColor = vec4(light.x + heat, light.y + (step*255 - heat/2), light.z + phiv, 1);
            // float dist_2 = phiv/4.0;
            gl_FragColor = vec4((light.x + heat*2.5) , (light.y - heat) + step, (light.z + phiv/3.0) - heat*2.5 - step, 1.0);
        }
        else{
            gl_FragColor = vec4((light.x + heat*2.5) , (light.y - heat), (light.z + phiv/3.0) - heat*2.5, 1.0);
        }
    }
    else{
        // vec3 dist = vec3(phiv, phiv, phiv);
        gl_FragColor = vec4((light.x + heat), light.y, light.z,  1.0);
    }
}
