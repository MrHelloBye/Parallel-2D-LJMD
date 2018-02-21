#version 400
in vec3 color;
out vec4 frag_color;
void main()
{
    //gl_FragColor = vec4(frag_color,1.0);
    frag_color = vec4(color,1.0);
}

