use std::ops::Deref;
use gl::types::*;
use std::str;
use std::fs;
use std::ffi::CString;
use std::ptr;

impl Deref for GLProgram {
    type Target = GLuint;
    fn deref(&self) -> &Self::Target {
        &self.handle
    }
}
pub struct GLProgram {
    handle : GLuint,
    vertex_shader_handle : GLuint,
    fragment_shader_handle : GLuint
}

impl Drop for GLProgram {
    fn drop(&mut self) {
        // Cleanup
        unsafe {
            gl::DeleteProgram(self.handle);
            gl::DeleteShader(self.vertex_shader_handle);
            gl::DeleteShader(self.fragment_shader_handle);
        }
    }
}

impl GLProgram {

    pub fn new(vertex_shader : &str, fragment_shader : &str) -> GLProgram {
        // Create GLSL shaders
        let vs = Self::compile_shader(&fs::read_to_string(vertex_shader).unwrap(), gl::VERTEX_SHADER);
        let fs = Self::compile_shader(&fs::read_to_string(fragment_shader).unwrap(), gl::FRAGMENT_SHADER);
        let program_handle = Self::link_program(vs, fs);

        GLProgram { handle : program_handle, fragment_shader_handle : fs, vertex_shader_handle : vs}
    }

fn compile_shader(src: &str, ty: GLenum) -> GLuint {
    let shader;
    unsafe {
        shader = gl::CreateShader(ty);
        // Attempt to compile the shader
        let c_str = CString::new(src.as_bytes()).unwrap();
        gl::ShaderSource(shader, 1, &c_str.as_ptr(), ptr::null());
        gl::CompileShader(shader);

        // Get the compile status
        let mut status = gl::FALSE as GLint;
        gl::GetShaderiv(shader, gl::COMPILE_STATUS, &mut status);

        // Fail on error
        if status != (gl::TRUE as GLint) {
            let mut len = 0;
            gl::GetShaderiv(shader, gl::INFO_LOG_LENGTH, &mut len);
            let mut buf = Vec::with_capacity(len as usize);
            buf.set_len((len as usize) - 1); // subtract 1 to skip the trailing null character
            gl::GetShaderInfoLog(
                shader,
                len,
                ptr::null_mut(),
                buf.as_mut_ptr() as *mut GLchar,
            );
            panic!(
                "{}",
                str::from_utf8(&buf)
                    .ok()
                    .expect("ShaderInfoLog not valid utf8")
            );
        }
    }
    shader
}

fn link_program(vs: GLuint, fs: GLuint) -> GLuint {
    unsafe {
        let program = gl::CreateProgram();
        gl::AttachShader(program, vs);
        gl::AttachShader(program, fs);
        gl::LinkProgram(program);
        // Get the link status
        let mut status = gl::FALSE as GLint;
        gl::GetProgramiv(program, gl::LINK_STATUS, &mut status);

        // Fail on error
        if status != (gl::TRUE as GLint) {
            let mut len: GLint = 0;
            gl::GetProgramiv(program, gl::INFO_LOG_LENGTH, &mut len);
            let mut buf = Vec::with_capacity(len as usize);
            buf.set_len((len as usize) - 1); // subtract 1 to skip the trailing null character
            gl::GetProgramInfoLog(
                program,
                len,
                ptr::null_mut(),
                buf.as_mut_ptr() as *mut GLchar,
            );
            panic!(
                "{}",
                str::from_utf8(&buf)
                    .ok()
                    .expect("ProgramInfoLog not valid utf8")
            );
        }
        program
    }
}
}