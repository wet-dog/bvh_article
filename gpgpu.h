#pragma once

namespace Tmpl8
{

// application class
class GPGPUApp : public TheApp
{
public:
	// game flow methods
	void Init();
	void AnimateScene();
	float3 Trace( Ray& ray, int rayDepth = 0 );
	void Tick( float deltaTime );
	void Shutdown() { /* implement if you want to do something on exit */ }
	// input handling
	void MouseUp( int button ) { /* implement if you want to detect mouse button presses */ }
	void MouseDown( int button ) { /* implement if you want to detect mouse button presses */ }
	void MouseMove( int x, int y ) { mousePos.x = x, mousePos.y = y; }
	void MouseWheel( float y ) { /* implement if you want to handle the mouse wheel */ }
	void KeyUp( int key ) { /* implement if you want to handle keys */ }
	void KeyDown( int key ) { /* implement if you want to handle keys */ }
	// data members
	int2 mousePos;
	Mesh* mesh;
	BVHInstance bvhInstance[256];
	TLAS tlas;
	float3* position, *direction, *orientation;
	float3 p0, p1, p2;	// virtual screen plane corners
	float* skyPixels;
	int skyWidth, skyHeight, skyBpp;
	Kernel* tracer;		// the ray tracing kernel
	Buffer* target;		// buffer encapsulating texture that holds the rendered image
	Buffer* skyData;	// buffer for the skydome texture
};

} // namespace Tmpl8