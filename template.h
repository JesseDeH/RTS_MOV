// Template, major revision 6, update for INFOMOV
// IGAD/NHTV/UU - Jacco Bikker - 2006-2016

#pragma once

#include "math.h"
#include "stdlib.h"
#include "emmintrin.h"
#include "stdio.h"
#include "windows.h"
#include "surface.h"

namespace Tmpl8 {

// vectors
class float2 // adapted from https://github.com/dcow/RayTracer
{
public:
	union { struct { float x, y; }; float cell[2]; };
	float2() {}
	float2( float v ) : x( v ), y( v ) {}
	float2( float x, float y ) : x( x ), y( y ) {}
	float2 operator - () const { return float2( -x, -y ); }
	float2 operator + ( const float2& addOperand ) const { return float2( x + addOperand.x, y + addOperand.y ); }
	float2 operator - ( const float2& operand ) const { return float2( x - operand.x, y - operand.y ); }
	float2 operator * ( const float2& operand ) const { return float2( x * operand.x, y * operand.y ); }
	float2 operator * ( float operand ) const { return float2( x * operand, y * operand ); }
	void operator -= ( const float2& a ) { x -= a.x; y -= a.y; }
	void operator += ( const float2& a ) { x += a.x; y += a.y; }
	void operator *= ( const float2& a ) { x *= a.x; y *= a.y; }
	void operator *= ( float a ) { x *= a; y *= a; }
	float& operator [] ( const int idx ) { return cell[idx]; }
	float length() { return sqrtf( x * x + y * y ); }
	float sqlength() { return x * x + y * y; }
	float2 normalized() { float r = 1.0f / length(); return float2( x * r, y * r ); }
	void normalize() { float r = 1.0f / length(); x *= r; y *= r; }
	static float2 normalize( float2 v ) { return v.normalized(); }
	float dot( const float2& operand ) const { return x * operand.x + y * operand.y; }
};

class float3
{
public:
	union { struct { float x, y, z; }; float cell[3]; };
	float3() {}
	float3( float v ) : x( v ), y( v ), z( v ) {}
	float3( float x, float y, float z ) : x( x ), y( y ), z( z ) {}
	float3 operator - () const { return float3( -x, -y, -z ); }
	float3 operator + ( const float3& addOperand ) const { return float3( x + addOperand.x, y + addOperand.y, z + addOperand.z ); }
	float3 operator - ( const float3& operand ) const { return float3( x - operand.x, y - operand.y, z - operand.z ); }
	float3 operator * ( const float3& operand ) const { return float3( x * operand.x, y * operand.y, z * operand.z ); }
	void operator -= ( const float3& a ) { x -= a.x; y -= a.y; z -= a.z; }
	void operator += ( const float3& a ) { x += a.x; y += a.y; z += a.z; }
	void operator *= ( const float3& a ) { x *= a.x; y *= a.y; z *= a.z; }
	void operator *= ( const float a ) { x *= a; y *= a; z *= a; }
	float operator [] ( const unsigned int& idx ) const { return cell[idx]; }
	float& operator [] ( const unsigned int& idx ) { return cell[idx]; }
	float length() const { return sqrtf( x * x + y * y + z * z ); }
	float sqrLentgh() const { return x * x + y * y + z * z; }
	float3 normalized() const { float r = 1.0f / length(); return float3( x * r, y * r, z * r ); }
	void normalize() { float r = 1.0f / length(); x *= r; y *= r; z *= r; }
	static float3 normalize( const float3 v ) { return v.normalized(); }
	float3 cross( const float3& operand ) const
	{
		return float3( y * operand.z - z * operand.y, z * operand.x - x * operand.z, x * operand.y - y * operand.x );
	}
	float dot( const float3& operand ) const { return x * operand.x + y * operand.y + z * operand.z; }
};

float3 normalize( const float3& v );
float2 normalize( const float2& v );
float3 cross( const float3& a, const float3& b );
float dot( const float3& a, const float3& b );
float dot( const float2& a, const float2& b );
float3 operator * ( const float& s, const float3& v );
// float2 operator * ( const float& s, const float2& v );
float3 operator * ( const float3& v, const float& s );
float2 operator * ( float2& v, float& s );
float length( const float2& v );
float sqlength(const float2& v);
unsigned int morton(const unsigned int x, const unsigned int y);
unsigned int splitBitwise(const unsigned int x);

};

#include "game.h"
#include <vector>
#include "freeimage.h"
#include "threads.h"

extern "C" 
{ 
#include "glew.h" 
}
#include "gl.h"
#include "io.h"
#include <ios>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "SDL.h"
#include "SDL_syswm.h"
#include "wglext.h"
#include "fcntl.h"

using namespace Tmpl8;				// to use template classes
using namespace std;				// to use stl vectors

inline float Rand( float range ) { return ((float)rand() / RAND_MAX) * range; }
inline int IRand( int range ) { return rand() % range; }
int filesize( FILE* f );

namespace Tmpl8 {

#define PI					3.14159265358979323846264338327950f
#define INVPI				0.31830988618379067153776752674503f

#define MALLOC64(x)			_aligned_malloc(x,64)
#define FREE64(x)			_aligned_free(x)
#define PREFETCH(x)			_mm_prefetch((const char*)(x),_MM_HINT_T0)
#define PREFETCH_ONCE(x)	_mm_prefetch((const char*)(x),_MM_HINT_NTA)
#define PREFETCH_WRITE(x)	_m_prefetchw((const char*)(x))
#define loadss(mem)			_mm_load_ss((const float*const)(mem))
#define broadcastps(ps)		_mm_shuffle_ps((ps),(ps), 0)
#define broadcastss(ss)		broadcastps(loadss((ss)))

struct Timer 
{ 
	typedef long long value_type; 
	static double inv_freq; 
	value_type start; 
	Timer() : start( get() ) { init(); } 
	float elapsed() const { return (float)((get() - start) * inv_freq); } 
	static value_type get() 
	{ 
		LARGE_INTEGER c; 
		QueryPerformanceCounter( &c ); 
		return c.QuadPart; 
	} 
	static double to_time(const value_type vt) { return double(vt) * inv_freq; } 
	void reset() { start = get(); }
	static void init() 
	{ 
		LARGE_INTEGER f; 
		QueryPerformanceFrequency( &f ); 
		inv_freq = 1000./double(f.QuadPart); 
	} 
}; 

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned char byte;

#define BADFLOAT(x) ((*(uint*)&x & 0x7f000000) == 0x7f000000)

}; // namespace Tmpl8