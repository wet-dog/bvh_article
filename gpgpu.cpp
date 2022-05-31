#include "precomp.h"
#include "bvh.h"
#include "gpgpu.h"

// THIS FILE IS UNDER CONSTRUCTION - WILL BE USED WITH ARTICLE #9

// THIS SOURCE FILE:
// Code for the article "How to Build a BVH", part 9: GPGPU.
// This version shows how to render a scene using ray tracing on the
// GPU - without hardware ray tracing. The scene (and accstruc) is
// fully maintained on the CPU.
// Feel free to copy this code to your own framework. Absolutely no
// rights are reserved. No responsibility is accepted either.
// For updates, follow me on twitter: @j_bikker.

TheApp* CreateApp() { return new GPGPUApp(); }

// GPGPUApp implementation

void GPGPUApp::Init()
{
	mesh = new Mesh( "assets/teapot.obj", "assets/bricks.png" );
	for (int i = 0; i < 16; i++)
		bvhInstance[i] = BVHInstance( mesh->bvh, i );
	tlas = TLAS( bvhInstance, 16 );
	// load HDR sky
	int bpp = 0;
	skyPixels = stbi_loadf( "assets/sky_19.hdr", &skyWidth, &skyHeight, &skyBpp, 0 );
	for (int i = 0; i < skyWidth * skyHeight * 3; i++) skyPixels[i] = sqrtf( skyPixels[i] );
	// prepare OpenCL
	tracer = new Kernel( "cl/kernels.cl", "render" );
	target = new Buffer( SCRWIDTH * SCRHEIGHT * 4 ); // intermediate screen buffer / render target
	skyData = new Buffer( skyWidth * skyHeight * 3 * sizeof( float ), skyPixels );
	skyData->CopyToDevice();
}

void GPGPUApp::AnimateScene()
{
	// animate the scene
	static float a[16] = { 0 }, h[16] = { 5, 4, 3, 2, 1, 5, 4, 3 }, s[16] = { 0 };
	for (int i = 0, x = 0; x < 4; x++) for (int y = 0; y < 4; y++, i++)
	{
		mat4 R, T = mat4::Translate( (x - 1.5f) * 2.5f, 0, (y - 1.5f) * 2.5f );
		if ((x + y) & 1) R = mat4::RotateY( a[i] );
		else R = mat4::Translate( 0, h[i / 2], 0 );
		if ((a[i] += (((i * 13) & 7) + 2) * 0.005f) > 2 * PI) a[i] -= 2 * PI;
		if ((s[i] -= 0.01f, h[i] += s[i]) < 0) s[i] = 0.2f;
		bvhInstance[i].SetTransform( T * R * mat4::Scale( 1.5f ) );
	}
	// update the TLAS
	tlas.Build();
}
  
void GPGPUApp::Tick( float deltaTime )
{
	// update the TLAS
	AnimateScene();
	// setup screen plane in world space
	static float angle = 0, ar = (float)SCRWIDTH / SCRHEIGHT; angle += 0.01f;
	mat4 M1 = mat4::RotateY( angle ), M2 = M1 * mat4::RotateX( -0.65f );
	p0 = TransformPosition( float3( -1 * ar, 1, 1.5f ), M2 );
	p1 = TransformPosition( float3( 1 * ar, 1, 1.5f ), M2 );
	p2 = TransformPosition( float3( -1 * ar, -1, 1.5f ), M2 );
	float3 camPos = TransformPosition( float3( 0, -2, -8.5f ), M1 );
	// render the scene using the GPU
	tracer->SetArguments( target, skyData, camPos, p0, p1, p2 );
	tracer->Run( SCRWIDTH * SCRHEIGHT );
	// obtain the rendered result
	target->CopyFromDevice();
	memcpy( screen->pixels, target->GetHostPtr(), target->size );
}

// EOF