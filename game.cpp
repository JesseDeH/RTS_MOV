#include "template.h"
#include <list>

// global data (source scope)
static Game* game;

// mountain peaks (push player away)
static float peakx[16] = { 248, 537, 695, 867, 887, 213, 376, 480, 683, 984, 364, 77,  85, 522, 414, 856 };
static float peaky[16] = { 199, 223, 83, 374,  694, 639, 469, 368, 545, 145,  63, 41, 392, 285, 447, 352 };
static float peakh[16] = { 200, 150, 160, 255, 200, 255, 200, 300, 120, 100,  80, 80,  80, 160, 160, 160 };

// player, bullet and smoke data
static int aliveP1 = MAXP1, aliveP2 = MAXP2;
static Bullet bullet[MAXBULLET];

int tankGrid[GRIDWIDTH * GRIDHEIGHT * GRIDROW];
float sinTable[720];
float cosTable[720];

// smoke particle effect tick function
void Smoke::Tick()
{
	unsigned int p = frame >> 3;
	if (frame < 64) if (!(frame++ & 7)) puff[p].x = xpos, puff[p].y = ypos << 8, puff[p].vy = -450, puff[p].life = 63;
	for ( unsigned int i = 0; i < p; i++ ) if ((frame < 64) || (i & 1))
	{
		puff[i].x++, puff[i].y += puff[i].vy, puff[i].vy += 3;
		int frame = (puff[i].life > 13)?(9 - (puff[i].life - 14) / 5):(puff[i].life / 2);
		game->m_Smoke->SetFrame( frame );
		game->m_Smoke->Draw( puff[i].x - 12, (puff[i].y >> 8) - 12, game->m_Surface );
		if (!--puff[i].life) puff[i].x = xpos, puff[i].y = ypos << 8, puff[i].vy = -450, puff[i].life = 63;
	}
}

// bullet Tick function
void Bullet::Tick()
{
	if (!(flags & Bullet::ACTIVE)) return;
	float2 prevpos = pos;
	pos += speed * 1.5f, prevpos -= pos - prevpos;
	game->m_Surface->AddLine( prevpos.x, prevpos.y, pos.x, pos.y, 0x555555 );
	if ((pos.x < 0) || (pos.x > (SCRWIDTH - 1)) || (pos.y < 0) || (pos.y > (SCRHEIGHT - 1))) flags = 0; // off-screen
	unsigned int start = 0, end = MAXP1;
	if (flags & P1) start = MAXP1, end = MAXP1 + MAXP2;
	for ( unsigned int i = start; i < end; i++ ) // check all opponents
	{
		Tank* t = game->m_Tank[i];
		if (!((t->flags & Tank::ACTIVE) && (pos.x > (t->pos.x - 2)) && (pos.y > (t->pos.y - 2)) && (pos.x < (t->pos.x + 2)) && (pos.y < (t->pos.y + 2)))) continue;
		if (t->flags & Tank::P1) aliveP1--; else aliveP2--; // update counters
		t->flags &= Tank::P1|Tank::P2;	// kill tank
		flags = 0;						// destroy bullet
		break;
	}
}

// Tank::Fire - spawns a bullet
void Tank::Fire( unsigned int party, float2& pos, float2& dir )
{
	for ( unsigned int i = 0; i < MAXBULLET; i++ ) if (!(bullet[i].flags & Bullet::ACTIVE))
	{
		bullet[i].flags |= Bullet::ACTIVE + party; // set owner, set active
		bullet[i].pos = pos, bullet[i].speed = speed;
		break;
	}
}

// Tank::Tick - update single tank
void Tank::Tick()
{
	if (!(flags & ACTIVE)) // dead tank
	{
		smoke.xpos = (int)pos.x, smoke.ypos = (int)pos.y;
		return smoke.Tick();
	}
	float2 force = normalize( target - pos );
	// evade mountain peaks
	for ( unsigned int i = 0; i < 16; i++ )
	{
		float2 d( pos.x - peakx[i], pos.y - peaky[i] );
		float sd = (d.x * d.x + d.y * d.y) * 0.2f;
		if (sd < 1500) 
		{
			force += d * 0.03f * (peakh[i] / sd);
			float r = sqrtf( sd );
			for( int j = 0; j < 720; j++ )
			{
#ifdef SINCOSLOOKUP
				float x = peakx[i] + r * sinTable[j];
				float y = peaky[i] + r * cosTable[j];
#else
				float x = peakx[i] + r * sinf( (float)j * PI / 360.0f );
				float y = peaky[i] + r * cosf( (float)j * PI / 360.0f );
#endif
				game->m_Surface->AddPlot( (int)x, (int)y, 0x000500 );
			}
		}
	}

	// evade other tanks
#ifdef GRID
	int grid_x = (int)(pos.x / GRIDSIZE);
	int grid_y = (int)(pos.y / GRIDSIZE);
	int gridPointer = ( grid_x + (grid_y * GRIDWIDTH) ) * GRIDROW;
	int gridCounter = tankGrid[gridPointer];
	int start_x = (grid_x - 1) < 0 ? grid_x : (grid_x - 1);
	int start_y = (grid_y - 1) < 0 ? grid_y : (grid_y - 1);
	int end_x = (grid_x + 1) >= GRIDWIDTH ? grid_x : (grid_x + 1);
	int end_y = (grid_y + 1) >= GRIDHEIGHT ? grid_y : (grid_y + 1);
	for (int i = start_x; i <= end_x; i++)
	{
		for (int j = start_y; j <= end_y; j++)
		{
			int c_x = i;
			int c_y = j;
			int c_gridPointer = (c_x + (c_y * GRIDWIDTH)) * GRIDROW;
			int c_gridCounter = tankGrid[c_gridPointer];
			for (int c = 1; c <= c_gridCounter; c++)
			{
				int c_tankPointer = tankGrid[c_gridPointer + c];
				if (game->m_Tank[c_tankPointer] == this) continue;
				float2 d = pos - game->m_Tank[c_tankPointer]->pos;
				if (length(d) < 8) force += normalize(d) * 2.0f;
				else if (length(d) < 16) force += normalize(d) * 0.4f;
			}
		}
	}


#else
	for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++)
	{
		if (game->m_Tank[i] == this) continue;
		float2 d = pos - game->m_Tank[i]->pos;
		if (length(d) < 8) force += normalize(d) * 2.0f;
		else if (length(d) < 16) force += normalize(d) * 0.4f;
	}
#endif // GRID

	// evade user dragged line
	if ((flags & P1) && (game->m_LButton))
	{
		float x1 = (float)game->m_DStartX, y1 = (float)game->m_DStartY;
		float x2 = (float)game->m_MouseX, y2 = (float)game->m_MouseY;
		float2 N = normalize( float2( y2 - y1, x1 - x2 ) );
		float dist = dot( N, pos ) - dot( N, float2( x1, y1 ) );
		if (fabs( dist ) < 10) if (dist > 0) force += N * 20; else force -= N * 20;
	}

	// update speed using accumulated force
	speed += force, speed = normalize( speed ), pos += speed * maxspeed * 0.5f;
#ifdef GRID
	UpdateGrid();
#endif // GRID

#ifndef GRID2
	// shoot, if reloading completed
	if (--reloading >= 0) return;
	unsigned int start = 0, end = MAXP1;
	if (flags & P1) start = MAXP1, end = MAXP1 + MAXP2;
	for (unsigned int i = start; i < end; i++) if (game->m_Tank[i]->flags & ACTIVE)
	{
		float2 d = game->m_Tank[i]->pos - pos;
		if ((length(d) < 100) && (dot(normalize(d), speed) > 0.99999f))
		{
			Fire(flags & (P1 | P2), pos, speed); // shoot
			reloading = 200; // and wait before next shot is ready
			break;
		}
	}
#endif // !GRID
}

void Tank::UpdateGrid()
{
	// Remove from old position
	int gridCounter = tankGrid[gridPointer];
	tankGrid[gridPosition] = tankGrid[gridPointer + gridCounter];
	tankGrid[gridPointer]--;
	
	// Add to new position
	ADDTOGRID();
}

void Tank::ADDTOGRID()
{
	int x = (int)(pos.x / GRIDSIZE);
	int y = (int)(pos.y / GRIDSIZE) * (GRIDWIDTH);
	gridPointer = (x+y) * GRIDROW;
	int gridCounter = ++tankGrid[gridPointer];
	if (gridCounter > GRIDROW)
		throw ERROR;
	gridPosition = gridPointer + gridCounter;
	tankGrid[gridPosition] = listPosition;
}

void Tank::CheckShooting()
{
#ifdef GRID2
	if (--reloading >= 0 || !(flags & ACTIVE)) return;
	int team = P1;
	if (flags & P1)
		team = P2;
	int grid_x = (int)(pos.x / GRIDSIZE);
	int grid_y = (int)(pos.y / GRIDSIZE);

	int start_x = max(grid_x - 7, 0);
	int start_y = max(grid_y - 7, 0);
	int end_x = min(grid_x + 7, GRIDWIDTH - 1);
	int end_y = min(grid_y + 7, GRIDHEIGHT - 1);

	for (int i = start_x; i <= end_x; i++)
	{
		for (int j = start_y; j <= end_y; j++)
		{
			int c_gridPointer = (i + (j * GRIDWIDTH)) * GRIDROW;
			int c_gridCounter = tankGrid[c_gridPointer];
			for (int c = c_gridPointer + 1; c <= c_gridPointer + c_gridCounter; c++)
			{
				int c_tankPointer = tankGrid[c];
				if (game->m_Tank[c_tankPointer]->flags & team)
				{
					if (game->m_Tank[c_tankPointer]->flags & ACTIVE)
					{
						float2 d = game->m_Tank[c_tankPointer]->pos - pos;
						if ((length(d) < 100) && (dot(normalize(d), speed) > 0.99999f))
						{
							Fire(flags & (P1 | P2), pos, speed); // shoot
							reloading = 200; // and wait before next shot is ready
							break;
						}
					}
				}
			}
		}
	}
#endif

#ifndef GRID2
	if (--reloading >= 0) return;
	unsigned int start = 0, end = MAXP1;
	if (flags & P1) start = MAXP1, end = MAXP1 + MAXP2;
	for (unsigned int i = start; i < end; i++) if (game->m_Tank[i]->flags & ACTIVE)
	{
		float2 d = game->m_Tank[i]->pos - pos;
		if ((length(d) < 100) && (dot(normalize(d), speed) > 0.99999f))
		{
			Fire(flags & (P1 | P2), pos, speed); // shoot
			reloading = 200; // and wait before next shot is ready
			break;
		}
	}
#endif
}

// Game::Init - Load data, setup playfield
void Game::Init()
{
	for (int i = 0; i < GRIDWIDTH*GRIDHEIGHT*GRIDROW; i++)
	{
		tankGrid[i] = 0;
	}
#ifdef SINCOSLOOKUP
	for (int j = 0; j < 720; j++)
	{
		sinTable[j] = sinf((float)j * PI / 360.0f);
		cosTable[j] = cosf((float)j * PI / 360.0f);
	}
#endif // SINCOSLOOKUP

	m_Heights = new Surface( "testdata/heightmap.png" ), m_Backdrop = new Surface( 1024, 768 ), m_Grid = new Surface( 1024, 768 );
	Pixel* a1 = m_Grid->GetBuffer(), *a2 = m_Backdrop->GetBuffer(), *a3 = m_Heights->GetBuffer();
	for ( int y = 0; y < 768; y++ ) for ( int idx = y * 1024, x = 0; x < 1024; x++, idx++ ) a1[idx] = (((x & 31) == 0) | ((y & 31) == 0)) ? 0x6600 : 0;
	for ( int y = 0; y < 767; y++ ) for ( int idx = y * 1024, x = 0; x < 1023; x++, idx++ ) 
	{
		float3 N = normalize( float3( (float)(a3[idx + 1] & 255) - (a3[idx] & 255), 1.5f, (float)(a3[idx + 1024] & 255) - (a3[idx] & 255) ) ), L( 1, 4, 2.5f );
		float h = (float)(a3[x + y * 1024] & 255) * 0.0005f, dx = x - 512.f, dy = y - 384.f, d = sqrtf( dx * dx + dy * dy ), dt = dot( N, normalize( L ) );
		int u = max( 0, min( 1023, (int)(x - dx * h) ) ), v = max( 0, min( 767, (int)(y - dy * h) ) ), r = (int)Rand( 255 );
		a2[idx] = AddBlend( a1[u + v * 1024], ScaleColor( ScaleColor( 0x33aa11, r ) + ScaleColor( 0xffff00, (255 - r) ), (int)(max( 0.0f, dt ) * 80.0f) + 10 ) );
	}
	m_Tank = new Tank*[MAXP1 + MAXP2];
	m_P1Sprite = new Sprite( new Surface( "testdata/p1tank.tga" ), 1, Sprite::FLARE );
	m_P2Sprite = new Sprite( new Surface( "testdata/p2tank.tga" ), 1, Sprite::FLARE );
	m_PXSprite = new Sprite( new Surface( "testdata/deadtank.tga" ), 1, Sprite::BLACKFLARE );
	m_Smoke = new Sprite( new Surface( "testdata/smoke.tga" ), 10, Sprite::FLARE );
	// create blue tanks
	for ( unsigned int i = 0; i < MAXP1; i++ )
	{
		Tank* t = m_Tank[i] = new Tank();
		t->pos = float2( (float)((i % 5) * 20), (float)((i / 5) * 20 + 50) );
#ifdef GRID
		t->listPosition = i;
		t->ADDTOGRID();
#endif //GRID
		t->target = float2( SCRWIDTH, SCRHEIGHT ); // initially move to bottom right corner
		t->speed = float2( 0, 0 ), t->flags = Tank::ACTIVE|Tank::P1, t->maxspeed = (i < (MAXP1 / 2))?0.65f:0.45f;
	}
	// create red tanks
	for ( unsigned int i = 0; i < MAXP2; i++ )
	{
		Tank* t = m_Tank[i + MAXP1] = new Tank();
		t->pos = float2( (float)((i % 12) * 20 + 900), (float)((i / 12) * 20 + 600) );
#ifdef GRID
		t->listPosition = i + MAXP1;
		t->ADDTOGRID();
#endif //GRID
		t->target = float2( 424, 336 ); // move to player base
		t->speed = float2( 0, 0 ), t->flags = Tank::ACTIVE|Tank::P2, t->maxspeed = 0.3f;
	}
	game = this; // for global reference
	m_LButton = m_PrevButton = false;
}

// Game::DrawTanks - draw the tanks
void Game::DrawTanks()
{
	for ( unsigned int i = 0; i < (MAXP1 + MAXP2); i++ )
	{
		Tank* t = m_Tank[i];
		float x = t->pos.x, y = t->pos.y;
		float2 p1( x + 70 * t->speed.x + 22 * t->speed.y, y + 70 * t->speed.y - 22 * t->speed.x );
		float2 p2( x + 70 * t->speed.x - 22 * t->speed.y, y + 70 * t->speed.y + 22 * t->speed.x );
		if (!(m_Tank[i]->flags & Tank::ACTIVE)) m_PXSprite->Draw( (int)x - 4, (int)y - 4, m_Surface ); // draw dead tank
		else if (t->flags & Tank::P1) // draw blue tank
		{
			m_P1Sprite->Draw( (int)x - 4, (int)y - 4, m_Surface );
			m_Surface->Line( x, y, x + 8 * t->speed.x, y + 8 * t->speed.y, 0x4444ff );
		}
		else // draw red tank
		{
			m_P2Sprite->Draw( (int)x - 4, (int)y - 4, m_Surface );
			m_Surface->Line( x, y, x + 8 * t->speed.x, y + 8 * t->speed.y, 0xff4444 );
		}
		if ((x >= 0) && (x < SCRWIDTH) && (y >= 0) && (y < SCRHEIGHT))
			m_Backdrop->GetBuffer()[(int)x + (int)y * SCRWIDTH] = SubBlend( m_Backdrop->GetBuffer()[(int)x + (int)y * SCRWIDTH], 0x030303 ); // tracks
	}
}

// Game::PlayerInput - handle player input
void Game::PlayerInput()
{
	if (m_LButton)
	{
		if (!m_PrevButton) m_DStartX = m_MouseX, m_DStartY = m_MouseY, m_DFrames = 0; // start line
		m_Surface->ThickLine( m_DStartX, m_DStartY, m_MouseX, m_MouseY, 0xffffff );
		m_DFrames++;
	}
	else
	{
		if ((m_PrevButton) && (m_DFrames < 15)) // new target location
			for ( unsigned int i = 0; i < MAXP1; i++ ) m_Tank[i]->target = float2( (float)m_MouseX, (float)m_MouseY );
		m_Surface->Line( 0, (float)m_MouseY, SCRWIDTH - 1, (float)m_MouseY, 0xffffff );
		m_Surface->Line( (float)m_MouseX, 0, (float)m_MouseX, SCRHEIGHT - 1, 0xffffff );
	}
	m_PrevButton = m_LButton;
}

// Game::Tick - main game loop
void Game::Tick( float a_DT )
{
#ifdef GRID
#endif // 
	POINT p;
	GetCursorPos( &p );
	ScreenToClient( FindWindow( NULL, "Template" ), &p );
	m_LButton = (GetAsyncKeyState( VK_LBUTTON ) != 0), m_MouseX = p.x, m_MouseY = p.y;
	m_Backdrop->CopyTo( m_Surface, 0, 0 );
	for ( unsigned int i = 0; i < (MAXP1 + MAXP2); i++ ) m_Tank[i]->Tick();
	for (unsigned int i = 0; i < (MAXP1 + MAXP2); i++) m_Tank[i]->CheckShooting();
	for ( unsigned int i = 0; i < MAXBULLET; i++ ) bullet[i].Tick();
	DrawTanks();
	PlayerInput();
	char buffer[128];
	if ((aliveP1 > 0) && (aliveP2 > 0))
	{
		sprintf( buffer, "blue army: %03i  red army: %03i", aliveP1, aliveP2 );
		return m_Surface->Print( buffer, 10, 10, 0xffff00 );
	}
	if (aliveP1 == 0) 
	{
		sprintf( buffer, "sad, you lose... red left: %i", aliveP2 );
		return m_Surface->Print( buffer, 200, 370, 0xffff00 );
	}
	sprintf( buffer, "nice, you win! blue left: %i", aliveP1 );
	m_Surface->Print( buffer, 200, 370, 0xffff00 );
}