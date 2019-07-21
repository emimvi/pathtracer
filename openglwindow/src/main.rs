extern crate gl;
extern crate glutin;
mod glprogram;

use gl::types::*;
use std::ffi::CString;
use std::ffi::CStr;
use std::mem;
use std::ptr;
use glprogram::*;

struct GLTexture {
    pub handle : GLuint, 
    width : usize,
    height : usize,
}
impl GLTexture {
    fn new() -> GLTexture {
        let mut handle = 0;
        let width = 2;
        let height = 2;

        unsafe {
            gl::ActiveTexture(gl::TEXTURE0);
            gl::GenTextures(1, &mut handle);
            gl::BindTexture(gl::TEXTURE_2D, handle);
            let t = _TEXTURE_TEST.as_ptr();
            gl::TexImage2D(gl::TEXTURE_2D, 0, gl::RGBA as i32, width as GLsizei, height as GLsizei, 0, gl::RGBA, gl::UNSIGNED_BYTE, t as *const _);

            //gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_WRAP_S, gl::CLAMP_TO_EDGE as i32);
            //gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_WRAP_T, gl::CLAMP_TO_EDGE as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MIN_FILTER, gl::NEAREST as i32);
            gl::TexParameteri(gl::TEXTURE_2D, gl::TEXTURE_MAG_FILTER, gl::NEAREST as i32);
        }

        GLTexture {
            handle,
            width,
            height
        }
    }
}

// Vertex data
static _TEXTURE_TEST : [GLuint; 4] = [0xFFFF0000, 0xFF00FF00, 0xFF0000FF, 0xFFFFFFFF];
static VERTEX_DATA: [GLfloat; 16] = [
 -0.5, -0.5,     0.0, 0.0,  
  0.5, -0.5,     1.0, 0.0,  
  0.5,  0.5,     1.0, 1.0,   
 -0.5,  0.5,     0.0, 1.0
 ];
static VERTEX_INDICIES: [GLuint; 6] = [0, 1, 2, 0, 2, 3];

struct GLMesh {
    vao : GLuint,
    vbo : GLuint,
    ebo : GLuint,
    num_triangles : usize
}

impl GLMesh {
    fn draw(&self) {
        unsafe {
            gl::BindVertexArray(self.vao);
            gl::DrawElements(gl::TRIANGLES, (self.num_triangles*3) as i32, gl::UNSIGNED_INT, ptr::null());
            gl::BindVertexArray(0);
        }
    }

    fn new(vertices : &[GLfloat], indicies : &[GLuint]) -> GLMesh {
        let mut vao = 0;
        let mut vbo = 0;
        let mut ebo = 0;
        debug_assert_eq!(indicies.len() % 3, 0);
        unsafe {
            // Create Vertex Array Object
            gl::GenVertexArrays(1, &mut vao);
            gl::BindVertexArray(vao);

            // Create a Vertex Buffer Object and copy the vertex data to it
            gl::GenBuffers(1, &mut vbo);
            gl::BindBuffer(gl::ARRAY_BUFFER, vbo);
            gl::BufferData(
                gl::ARRAY_BUFFER,
                (vertices.len() * mem::size_of::<GLfloat>()) as GLsizeiptr,
                mem::transmute(&vertices[0]),
                gl::STATIC_DRAW,
            );

            gl::GenBuffers(1, &mut ebo);
            gl::BindBuffer(gl::ELEMENT_ARRAY_BUFFER, ebo);
            gl::BufferData(
                gl::ELEMENT_ARRAY_BUFFER,
                (indicies.len() * mem::size_of::<GLuint>()) as GLsizeiptr,
                mem::transmute(&indicies[0]),
                gl::STATIC_DRAW,
            );


            gl::EnableVertexAttribArray(0);
            gl::VertexAttribPointer(0, 2, gl::FLOAT, gl::FALSE as GLboolean, (mem::size_of::<GLfloat>()*4) as i32, ptr::null());
            gl::EnableVertexAttribArray(1);
            gl::VertexAttribPointer(1, 2, gl::FLOAT, gl::FALSE as GLboolean, (mem::size_of::<GLfloat>()*4) as i32, (mem::size_of::<GLfloat>()*2) as *const _);
            gl::BindBuffer(gl::ARRAY_BUFFER, 0);
            gl::BindVertexArray(0);
        }

        GLMesh {
            vao,
            vbo,
            ebo,
            num_triangles : indicies.len()/3
        }
    }
}

impl Drop for GLMesh {
    fn drop(&mut self) {
        //Cleanup
        unsafe {
                gl::DeleteBuffers(2, [self.vbo, self.ebo].as_ptr());
                gl::DeleteVertexArrays(1, &self.vao);
        }
    }
}

fn main() {
    let mut events_loop = glutin::EventsLoop::new();
    let window = glutin::WindowBuilder::new().with_dimensions(glutin::dpi::LogicalSize::new(512.0, 512.0));
    let context = glutin::ContextBuilder::new().build_windowed(window, &events_loop).unwrap();

    // It is essential to make the context current before calling `gl::load_with`.
    let gl_window = unsafe { context.make_current() }.unwrap();

    // Load the OpenGL function pointers
    // TODO: `as *const _` will not be needed once glutin is updated to the latest gl version
    gl::load_with(|symbol| gl_window.get_proc_address(symbol) as *const _);

    let program = GLProgram::new("src/vert.vert", "src/frag.frag");

    print_gl_info();

    let mesh;
    let tex : GLTexture;
    unsafe {
        gl::DebugMessageCallback(debug_output_gl, ptr::null());
        gl::Enable(gl::DEBUG_OUTPUT_SYNCHRONOUS);
        
        gl::UseProgram(*program);
        mesh = GLMesh::new(&VERTEX_DATA, &VERTEX_INDICIES);
        // Use shader program
        gl::BindFragDataLocation(*program, 0, CString::new("out_color").unwrap().as_ptr());

        // Specify the layout of the vertex data
        //let pos_attr = gl::GetAttribLocation(*program, CString::new("position").unwrap().as_ptr());

        tex = GLTexture::new();
        //let tex_loc = gl::GetUniformLocation(*program, CString::new("te").unwrap().as_ptr());
        //gl::Uniform1i(tex_loc, tex.handle as GLint);
        //println!("{}", tex_loc);
    }
    events_loop.run_forever(|event| {
        use glutin::{ControlFlow, Event, WindowEvent};

        if let Event::WindowEvent { event, .. } = event {
            if let WindowEvent::CloseRequested = event {
                return ControlFlow::Break;
            }
        }

        unsafe {
            // Clear the screen to black
            gl::ClearColor(0.3, 0.3, 0.3, 1.0);
            gl::Clear(gl::COLOR_BUFFER_BIT);

            //gl::ActiveTexture(gl::TEXTURE0);
            //gl::BindTexture(gl::TEXTURE_2D, tex.handle);
            mesh.draw();
        }

        gl_window.swap_buffers().unwrap();

        ControlFlow::Continue
    });

}

extern "system" fn debug_output_gl(_source: GLenum, _type: GLenum, _id: GLuint, _severity: GLenum, _length: GLsizei, message: *const GLchar, _user_param: *mut GLvoid) {
    println!("GL error: {:?}", unsafe { CString::from_raw(message as *mut _) });
}

fn print_gl_info() {
    unsafe {
    let what = CStr::from_ptr(gl::GetString(gl::VENDOR) as *const i8);
    let what1 = CStr::from_ptr(gl::GetString(gl::RENDERER) as *const i8);
    let what2 = CStr::from_ptr(gl::GetString(gl::VERSION) as *const i8);
    let what3 = CStr::from_ptr(gl::GetString(gl::SHADING_LANGUAGE_VERSION) as *const i8);

    println!("{:?}", what);
    println!("{:?}", what1);
    println!("{:?}", what2);
    println!("{:?}", what3);
    }
}