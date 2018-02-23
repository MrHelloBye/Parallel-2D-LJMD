#version 110
attribute vec2 vertex_position;
attribute vec2 circle_position;
attribute vec3 circle_color;
varying vec3 color;
void main()
{
    gl_Position = vec4(vertex_position+circle_position,0,1.0);
    color = circle_color;
}

