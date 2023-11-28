#include "precomp.h"
#include "basics.h"

// THIS SOURCE FILE:
// Code for the article "How to Build a BVH", part 1: basics. Link:
// https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics
// This is bare-bones BVH construction and traversal code, running in
// a minimalistic framework.
// Feel free to copy this code to your own framework. Absolutely no
// rights are reserved. No responsibility is accepted either.
// For updates, follow me on twitter: @j_bikker.

namespace
{
	struct TriangleIndices
	{
		std::vector<uint32_t> vertices = std::vector<uint32_t>(3);
		std::vector<uint32_t> uvs = std::vector<uint32_t>(3);
		std::vector<uint32_t> normals = std::vector<uint32_t>(3);
	};

	void DecrementIndex(std::vector<uint32_t>& indices)
	{
		for (auto& i : indices)
		{
			--i;
			assert(i >= 0);
		}
	}

	void ZeroIndex(TriangleIndices& triangleIndices)
	{
		DecrementIndex(triangleIndices.vertices);
		DecrementIndex(triangleIndices.uvs);
		DecrementIndex(triangleIndices.normals);
	};
}

TheApp* CreateApp() { return new BasicBVHApp(); }

// triangle count
#define COUNT	240

// forward declarations
void Subdivide( uint nodeIdx );
void UpdateNodeBounds( uint nodeIdx );

// minimal structs
struct Tri { float3 vertex0, vertex1, vertex2; float3 centroid; };
__declspec(align(32)) struct BVHNode
{
	float3 aabbMin, aabbMax;
	uint leftFirst, triCount;
	bool isLeaf() { return triCount > 0; }
};
struct Ray { float3 O, D; float t = 1e30f; };

// application data
Tri tri[COUNT];
uint triIdx[COUNT];
BVHNode bvhNode[COUNT * 2];
uint rootNodeIdx = 0, nodesUsed = 1;
int triCount = 0;

// functions

void IntersectTri( Ray& ray, const Tri& tri )
{
	const float3 edge1 = tri.vertex1 - tri.vertex0;
	const float3 edge2 = tri.vertex2 - tri.vertex0;
	const float3 h = cross( ray.D, edge2 );
	const float a = dot( edge1, h );
	if (a > -0.0001f && a < 0.0001f) return; // ray parallel to triangle
	const float f = 1 / a;
	const float3 s = ray.O - tri.vertex0;
	const float u = f * dot( s, h );
	if (u < 0 || u > 1) return;
	const float3 q = cross( s, edge1 );
	const float v = f * dot( ray.D, q );
	if (v < 0 || u + v > 1) return;
	const float t = f * dot( edge2, q );
	if (t > 0.0001f) ray.t = min( ray.t, t );
}

bool IntersectAABB( const Ray& ray, const float3 bmin, const float3 bmax )
{
	float tx1 = (bmin.x - ray.O.x) / ray.D.x, tx2 = (bmax.x - ray.O.x) / ray.D.x;
	float tmin = min( tx1, tx2 ), tmax = max( tx1, tx2 );
	float ty1 = (bmin.y - ray.O.y) / ray.D.y, ty2 = (bmax.y - ray.O.y) / ray.D.y;
	tmin = max( tmin, min( ty1, ty2 ) ), tmax = min( tmax, max( ty1, ty2 ) );
	float tz1 = (bmin.z - ray.O.z) / ray.D.z, tz2 = (bmax.z - ray.O.z) / ray.D.z;
	tmin = max( tmin, min( tz1, tz2 ) ), tmax = min( tmax, max( tz1, tz2 ) );
	return tmax >= tmin && tmin < ray.t && tmax > 0;
}

void IntersectBVH( Ray& ray, const uint nodeIdx )
{
	BVHNode& node = bvhNode[nodeIdx];
	if (!IntersectAABB( ray, node.aabbMin, node.aabbMax )) return;
	if (node.isLeaf())
	{
		for (uint i = 0; i < node.triCount; i++ )
			IntersectTri( ray, tri[triIdx[node.leftFirst + i]] );
	}
	else
	{
		IntersectBVH( ray, node.leftFirst );
		IntersectBVH( ray, node.leftFirst + 1 );
	}
}

void BuildBVH()
{
	// populate triangle index array
	for (int i = 0; i < COUNT; i++) triIdx[i] = i;
	// calculate triangle centroids for partitioning
	for (int i = 0; i < COUNT; i++)
		tri[i].centroid = (tri[i].vertex0 + tri[i].vertex1 + tri[i].vertex2) * 0.3333f;
	// assign all triangles to root node
	BVHNode& root = bvhNode[rootNodeIdx];
	root.leftFirst = 0, root.triCount = COUNT;
	UpdateNodeBounds( rootNodeIdx );
	// subdivide recursively
	Subdivide( rootNodeIdx );
}

void UpdateNodeBounds( uint nodeIdx )
{
	BVHNode& node = bvhNode[nodeIdx];
	node.aabbMin = float3( 1e30f );
	node.aabbMax = float3( -1e30f );
	for (uint first = node.leftFirst, i = 0; i < node.triCount; i++)
	{
		uint leafTriIdx = triIdx[first + i];
		Tri& leafTri = tri[leafTriIdx];
		node.aabbMin = fminf( node.aabbMin, leafTri.vertex0 ),
		node.aabbMin = fminf( node.aabbMin, leafTri.vertex1 ),
		node.aabbMin = fminf( node.aabbMin, leafTri.vertex2 ),
		node.aabbMax = fmaxf( node.aabbMax, leafTri.vertex0 ),
		node.aabbMax = fmaxf( node.aabbMax, leafTri.vertex1 ),
		node.aabbMax = fmaxf( node.aabbMax, leafTri.vertex2 );
	}
}

void Subdivide( uint nodeIdx )
{
	// terminate recursion
	BVHNode& node = bvhNode[nodeIdx];
	if (node.triCount <= 2) return;
	// determine split axis and position
	float3 extent = node.aabbMax - node.aabbMin;
	int axis = 0;
	if (extent.y > extent.x) axis = 1;
	if (extent.z > extent[axis]) axis = 2;
	float splitPos = node.aabbMin[axis] + extent[axis] * 0.5f;
	// in-place partition
	int i = node.leftFirst;
	int j = i + node.triCount - 1;
	while (i <= j)
	{
		if (tri[triIdx[i]].centroid[axis] < splitPos)
			i++;
		else
			swap( triIdx[i], triIdx[j--] );
	}
	// abort split if one of the sides is empty
	int leftCount = i - node.leftFirst;
	if (leftCount == 0 || leftCount == node.triCount) return;
	// create child nodes
	int leftChildIdx = nodesUsed++;
	int rightChildIdx = nodesUsed++;
	bvhNode[leftChildIdx].leftFirst = node.leftFirst;
	bvhNode[leftChildIdx].triCount = leftCount;
	bvhNode[rightChildIdx].leftFirst = i;
	bvhNode[rightChildIdx].triCount = node.triCount - leftCount;
	node.leftFirst = leftChildIdx;
	node.triCount = 0;
	UpdateNodeBounds( leftChildIdx );
	UpdateNodeBounds( rightChildIdx );
	// recurse
	Subdivide( leftChildIdx );
	Subdivide( rightChildIdx );
}

void BasicBVHApp::Init()
{
	std::string filename = "assets/teapot-low.obj";
	std::ifstream file(filename);
	if (file.fail()) return;

	std::array<float2, COUNT> UV = {}; // enough for dragon.obj
	std::array<float3, COUNT> N = {};
	std::array<float3, COUNT> P = {};
	int UVs = 0, Ns = 0, Ps = 0, a, b, c, d, e, f, g, h, i;
	for (std::string line; std::getline(file, line);)
	{
		const std::string header = line.substr(0, line.find(" "));
		if (header == "v")
		{
			const int ret = std::sscanf(line.c_str(), "v %f %f %f\n", &P[Ps].x, &P[Ps].y, &P[Ps].z);
			assert(ret == 3);
			++Ps;
		}

		if (header != "f")
		{
			continue;
		}
		else
		{
			TriangleIndices tri1;
			TriangleIndices tri2;

			// Counter-clockwise order
			const int ret = std::sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n",
				&tri1.vertices[0], &tri1.uvs[0], &tri1.normals[0],
				&tri1.vertices[1], &tri1.uvs[1], &tri1.normals[1],
				&tri1.vertices[2], &tri1.uvs[2], &tri1.normals[2],
				&tri2.vertices[2], &tri2.uvs[2], &tri2.normals[2]);
			assert(ret == 9 || ret == 12);

			ZeroIndex(tri1);
			ZeroIndex(tri2);

			if (ret == 9)
			{
				tri[triCount].vertex0 = P[tri1.vertices[0]];
				tri[triCount].vertex1 = P[tri1.vertices[1]];
				tri[triCount++].vertex2 = P[tri1.vertices[2]];
			}
			else if (ret == 12)
			{
				tri[triCount].vertex0 = P[tri1.vertices[0]];
				tri[triCount].vertex1 = P[tri1.vertices[1]];
				tri[triCount++].vertex2 = P[tri1.vertices[2]];
				tri[triCount].vertex0 = P[tri1.vertices[0]];
				tri[triCount].vertex1 = P[tri1.vertices[2]];
				tri[triCount++].vertex2 = P[tri2.vertices[2]];
			}
		}
	}

	// construct the BVH
	BuildBVH();
}

void BasicBVHApp::Tick( float deltaTime )
{
	// draw the scene
	screen->Clear( 0 );
	// define the corners of the screen in worldspace
	float3 p0( -1, 1, -15 ), p1( 1, 1, -15 ), p2( -1, -1, -15 );
	Ray ray;
	Timer t;
	for (int y = 0; y < SCRHEIGHT; y++) for (int x = 0; x < SCRWIDTH; x++)
	{
		// calculate the position of a pixel on the screen in worldspace
		float3 pixelPos = p0 + (p1 - p0) * (x / (float)SCRWIDTH) + (p2 - p0) * (y / (float)SCRHEIGHT);
		// define the ray in worldspace
		ray.O = float3( 0, 0, -18 );
		ray.D = normalize( pixelPos - ray.O );
		// initially the ray has an 'infinite length'
		ray.t = 1e30f;
	#if 0
		for( int i = 0; i < N; i++ ) IntersectTri( ray, tri[i] );
	#else
		IntersectBVH( ray, rootNodeIdx );
	#endif
		if (ray.t < 1e30f) screen->Plot( x, y, 0xffffff );
	}
	float elapsed = t.elapsed() * 1000;
	printf( "tracing time: %.2fms (%5.2fK rays/s)\n", elapsed, sqr( 630 ) / elapsed );
}

// EOF