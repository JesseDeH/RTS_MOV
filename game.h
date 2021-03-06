#ifndef I_GAME_H
#define I_GAME_H

#define	SCRWIDTH	1024
#define SCRHEIGHT	768

namespace Tmpl8 {

#define MAXP1		200	// increase to test your optimized code
#define MAXP2		 (4 * MAXP1)	// because the player is smarter than the AI
#define MAXBULLET	200

#define GRID
#define SINCOSLOOKUP
#define POSINGRID
	//#define MORTON
#define MOUNTAIN

#define GRIDSIZE		16												// How many pixels is one grid space
#define GRIDWIDTH		((11 * 20 + 900) / GRIDSIZE	+ 1)				// How many grid spaces wide is the grid
#define GRIDHEIGHT		((int) ((MAXP2 / 12) * 20 + 600) / GRIDSIZE + 1)	// How many grid spaces high is the grid
#ifdef POSINGRID
	#define GRIDROW			(32 * 3)												// How much data is there for era
#else
	#define GRIDROW			32
#endif

class Smoke
{
public:
	struct Puff { int x, y, vy, life; };
	Smoke() : active( false ), frame( 0 ) {};
	void Tick();
	Puff puff[8];
	bool active;
	int frame, xpos, ypos;
};

class Tank
{
public:
	enum { ACTIVE = 1, P1 = 2, P2 = 4 };
	Tank() : pos( float2( 0, 0 ) ), speed( float2( 0, 0 ) ), target( float2( 0, 0 ) ), reloading( 0 ) {};
	~Tank();
	void Fire( unsigned int party, float2& pos, float2& dir );
	void CheckShooting();
	void Tick();
	void UpdateGrid();
	void ADDTOGRID();
	float2 pos, speed, target;
	float maxspeed;
	int flags, reloading;
	int listPosition; //position in m_Tank array
	int gridPosition; //position of gridCounter in grid
	int gridPointer;  //position in grid
	Smoke smoke;
};

class Bullet
{
public:
	enum { ACTIVE = 1, P1 = 2, P2 = 4 };
	Bullet() : flags( 0 ) {};
	void Tick();
	float2 pos, speed;
	int flags;
};

class Surface;
class Surface8;
class Sprite;
class Game
{
public:
	void SetTarget( Surface* a_Surface ) { m_Surface = a_Surface; }
	void MouseMove( int x, int y ) { m_MouseX = x; m_MouseY = y; }
	void MouseButton( bool b ) { m_LButton = b; }
	void Init();
	void UpdateTanks();
	void UpdateBullets();
	void DrawTanks();
	void PlayerInput();
	void Tick( float a_DT );
	Surface* m_Surface, *m_Backdrop, *m_Heights, *m_Grid;
	unsigned int Morton(unsigned int x, unsigned int y);
	unsigned int intersperse(unsigned int x);
	Sprite* m_P1Sprite, *m_P2Sprite, *m_PXSprite, *m_Smoke;
	int m_ActiveP1, m_ActiveP2;
	int m_MouseX, m_MouseY, m_DStartX, m_DStartY, m_DFrames;
	bool m_LButton, m_PrevButton;
	Tank** m_Tank;
};

}; // namespace Templ8

#endif