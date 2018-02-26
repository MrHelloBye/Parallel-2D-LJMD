#version 400
in vec2 vertex_position;
in vec2 circle_position;
in vec3 circle_color;
out vec3 color;
void main()
{
    gl_Position = vec4(vertex_position+circle_position,0,1.0);
    color = circle_color;

}

