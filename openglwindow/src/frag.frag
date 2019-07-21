#version 450
out vec4 out_color;

in vec2 tex_coord;

uniform sampler2D tex;

void main() {
    out_color = texture(tex, tex_coord);
}