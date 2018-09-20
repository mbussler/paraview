#include "TimeStepControllerCPU.h"

//#include <boost/thread.hpp>
//#include <boost/bind.hpp>

#include "time.h"

TimeStepControllerCPU::TimeStepControllerCPU( )
{
	time=0.0;
	timeMax=0.0;
	currentTimestep=0;
	nextTimestep=0;
	overnextTimestep=0;
	integration=Dopri5_ATS;
	timestepCount=0;

	m_ats_params.tol = 10e-4f;
	m_ats_params.h_min = 0.001f;
	m_ats_params.rho = 0.9f;
	m_ats_params.eta = 2.0f;
};

TimeStepControllerCPU::~TimeStepControllerCPU()
{
	timesteps.clear();
};

void TimeStepControllerCPU::addTimestep( TimestepCPUPtr ts, int& timestepIndex )
{
	timesteps.push_back(ts);
	timestepIndex = currentTimestep;

	currentTimestep++;				
	timestepCount++;
};

/// Register a Timestep for a certain point-in-time
void TimeStepControllerCPU::registerTimestep( int ts_index, double ts_time  )
{
	assert( ts_index < timestepCount );
	
	TimePoint tp;
	tp.index = ts_index;
	tp.time = ts_time;
	timepoints.push_back(tp);
	timepoints.sort(TimePointSortPredicate);
	
	if( ts_time > timeMax )
		timeMax = ts_time;
};

TimestepCPUPtr TimeStepControllerCPU::getTimestep( int timestepIndex )
{
	assert( timestepIndex >=0 &&  timestepIndex < timestepCount);
	return timesteps[timestepIndex];
};

TimestepCPUPtr TimeStepControllerCPU::getCurrentTimestep()
{
	return timesteps[currentTimestep];
};

TetMeshGLPtr TimeStepControllerCPU::getCurrentMeshGL()
{
	return timesteps[currentTimestep]->getMeshGL();
};

TetMeshGLPtr TimeStepControllerCPU::getNextMeshGL()
{
	return timesteps[nextTimestep]->getMeshGL();
};

TimestepCPUPtr TimeStepControllerCPU::getNextTimestep()
{
	return timesteps[nextTimestep];
};
TimestepCPUPtr TimeStepControllerCPU::getOvernextTimestep()
{
	return timesteps[overnextTimestep];
};

void TimeStepControllerCPU::initTimesteps()
{
	for( int i=0; i< timesteps.size(); i++)
	{
		timesteps[i]->init(i);
	}
	
	resetTimepoints();

	getCurrentTimestep()->storeTraversedCells(true);
};

void TimeStepControllerCPU::resetTime()
{
	// reset time to the beginning
	time = 0.0f;

	// reset timepoints
	resetTimepoints();
};
void TimeStepControllerCPU::resetTimepoints()
{
	// set indices for current, next and overnext timestep
	currentTp = timepoints.begin();
	currentTimestep = currentTp->index;

	nextTp = currentTp;
	safeIncrement(timepoints, nextTp);
	nextTimestep = nextTp->index;

	overnextTp = nextTp;
	safeIncrement(timepoints, overnextTp);
	overnextTimestep = overnextTp->index;

};

void TimeStepControllerCPU::incrementTimepoints()
{
	// increment all indices
	currentTp = nextTp;
	currentTimestep = nextTimestep;

	nextTp = overnextTp;
	nextTimestep = overnextTimestep;

	safeIncrement( timepoints, overnextTp );
	overnextTimestep = overnextTp->index;

};

//	// store the traversed cells of the overnext ts for the kd-walk vis
//	timesteps[overnextTimestep]->storeTraversedCells( true );

void TimeStepControllerCPU::advectParticles( ParticleList& particles, float h )
{

	if( time+h < nextTp->time )
	{
		// calculate advection for current and next timestep
		calculateAdvection( particles, h, currentTimestep, nextTimestep );

		time += h;
		return;
	}
	else
	{
		// calculate Advection for the first half of the timestep
		float h1 = nextTp->time - time;

		calculateAdvection( particles, h1, currentTimestep, nextTimestep );
		time += h1;

		float h2 = h-h1;

		if( time < timeMax )
		{
			// increment timePoints
			incrementTimepoints();

			// locate all particles in the new timestep
			//getNextTimestep()->getCudaMesh()->resetTraversedCells();
			getNextTimestep()->locateParticles( particles );

			if( h2 > 0.05 * h )
				calculateAdvection( particles, h2, currentTimestep, nextTimestep );
		}
		else
		{
			// time has crossed the maximum time -> reset everything
			time = time - timeMax;

			// reset timepoints
			resetTimepoints();

			// locate particles
			getCurrentTimestep()->locateParticles( particles );
			getNextTimestep()->locateParticles( particles );

			if( h2 > 0.05 * h )
				calculateAdvection( particles, h2, currentTimestep, nextTimestep );

		}
		time += h2;
	}
}

void TimeStepControllerCPU::calculateAdvection( ParticleList& particles, float timestep, int ts_a, int ts_b )
{
	// The normalized time for the interpolation
	double a = (time - currentTp->time) / (nextTp->time - currentTp->time);

	// The normalized timestep
	double ah = timestep / (nextTp->time - currentTp->time);

	// Clear traversed cells
	timesteps[ts_a]->getMeshGL()->resetTraversedCells();
	timesteps[ts_b]->getMeshGL()->resetTraversedCells();

  timer1.restart();
	for( ParticleListIterator p = particles.begin(); p != particles.end(); p++)
	{
		// calculate velocity at particles position
		Vector velocity;

		CalculateIntegration(*p, velocity, timestep, a, ah, ts_a, ts_b );

		p->move( velocity * timestep );

		p->bOutOfField = ( p->Cell[currentTimestep] == -1 || p->Cell[nextTimestep] == -1);
	}
  timer1.stop();
}

void TimeStepControllerCPU::switchIntegration()
{
	switch( integration)
	{
	case Euler:
		integration = RK3;
		//printf("RK-3:");
		break;
	case RK3:
		integration = RK4;
		//printf("RK-4:");
		break;
	case RK4:
		integration = Dopri5;
		//printf("Dopri-5:");
		break;
	case Dopri5:
		integration = Dopri5_ATS;
		//printf("Dopri-5 ATS:");
		break;
	case Dopri5_ATS:
		integration = Euler;
		//printf("Euler:");
		break;
	}
};

void TimeStepControllerCPU::CalculateIntegration(Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b)
{
	Vector v_rk4;

	switch( integration )
	{
		default:
		case Euler:
			CalculateEulerIntegration( p, v, h, a, ah, ts_a, ts_b);
			break;
		case RK3:
			CalculateRK3Integration( p, v, h, a, ah, ts_a, ts_b);
			break;
		case RK4:
			CalculateRK4Integration( p, v, h, a, ah, ts_a, ts_b);
			break;
		case Dopri5:
			CalculateDopri5Integration( p, v_rk4, v, h, a, ah, ts_a, ts_b);
			break;
		case Dopri5_ATS:
			CalculateDopri5Integration_ATS( p, v, h, a, ah, ts_a, ts_b);
			break;
	}
};

void TimeStepControllerCPU::CalculateEulerIntegration(Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b)
{
	Vector v1, v2;
	timesteps[ts_a]->interpolateVelocityAt(p.Position, p.Cell[ts_a], v1);
	timesteps[ts_b]->interpolateVelocityAt(p.Position, p.Cell[ts_b], v2);

	v = (1-a)*v1 + a*v2;
};

void TimeStepControllerCPU::CalculateRK3Integration(Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b)
{
	Vector kx1, kx2;
	Vector y0, yx;

	Vector k1,k2,k3;

	const float h_div_2 = h/2.0;

	y0 = Pt2Vec(p.Position);

		timesteps[ts_a]->interpolateVelocityAt( y0.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( y0.pt(), p.Cell[ts_b], kx2 );

	float ax = a;
	k1 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h_div_2*k1;

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah/2.0f;
	k2 = (1-ax)*kx1 + ax*kx2;

	yx = y0 - h*k1 + 2.0*h*k2;

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah;
	k3 = (1-ax)*kx1 + ax*kx2;
		
	v = 1.0/6.0*k1  + 4.0/6.0*k2 + 1.0/6.0*k3;
};

void TimeStepControllerCPU::CalculateRK4Integration(Particle& p, Vector& v, float h, float a, float ah, int ts_a, int ts_b)
{
	Vector kx1, kx2;
	Vector y0, yx;

	Vector k1,k2,k3,k4;

	const float h_div_2 = h/2.0;

	y0 = Pt2Vec(p.Position);

		timesteps[ts_a]->interpolateVelocityAt( y0.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( y0.pt(), p.Cell[ts_b], kx2 );

	// c_j=0
	float ax = a;
	k1 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h_div_2*k1;

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );
	
	// c_j=1/2
	ax = a + ah/2.0;
	k2 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h_div_2*k2;

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	// c_j=1/2
	ax = a + ah/2.0;
	k3 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*k3;

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	// c_j=1
	ax = a + ah;
	k4 = (1-ax)*kx1 + ax*kx2;

	v = 1.0/6.0*k1 + 1.0/3.0*k2 + 1.0/3.0*k3 + 1.0/6.0*k4;
};

void TimeStepControllerCPU::CalculateDopri5Integration( Particle& p, Vector& v_rk4, Vector& v_rk5, float h, float a, float ah, int ts_a, int ts_b)
{
	float ax;
	Vector y0;
	Vector k0, k1, k2, k3, k4, k5, k6;
	Vector yx, kx1, kx2;

	// flow field evaluation in mesh 1 and mesh 2
	y0 = Pt2Vec(p.Position);

		timesteps[ts_a]->interpolateVelocityAt( y0.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( y0.pt(), p.Cell[ts_b], kx2 );

	// temporal interpolation
	ax = a;
	k0 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*(k0/5.0f);

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah/5.0f;
	k1 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*( 	3.0f/40.0f*k0 +
					9.0f/40.0f*k1 );

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah*(3.0f/10.0f );
	k2 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*(	44.0f/45.0f * k0 -
					56.0f/15.0f * k1 +
					32.0f/9.0f  * k2 );

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah*(4.0f/5.0f);
	k3 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*(	19372.0f/6561.0f  * k0 -
					25360.0f/2187.0f  * k1 +
					64448.0f/6561.0f  * k2 -
					212.0f/729.0f     * k3 );

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah*(8.0f/ 9.0f);
	k4 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*(	9017.0f/3168.0f   * k0 -
					355.0f/33.0f      * k1 +
					46732.0f/5247.0f  * k2 +
					49.0f/176.0f      * k3 -
					5103.0f/18656.0f  * k4);

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah;
	k5 = (1-ax)*kx1 + ax*kx2;

	yx = y0 + h*(	35.0f/384.0f     * k0 +
					//0.0f           * k1 +
					500.0f/1113.0f   * k2 +
					125.0f/192.0f    * k3 -
					2187.0f/6784.0f  * k4 +
					11.0f/84.0f      * k5 );

		timesteps[ts_a]->interpolateVelocityAt( yx.pt(), p.Cell[ts_a], kx1 );
		timesteps[ts_b]->interpolateVelocityAt( yx.pt(), p.Cell[ts_b], kx2 );

	ax = a + ah;
	k6 = (1-ax)*kx1 + ax*kx2;

	v_rk4 = \
		5179.0f / 57600.0f 	* k0 +
		// 0.0f 			* k1 +
		7571.0f / 16695.0f 	* k2 +
		393.0f / 640.0f 	* k3 -
		92097.0f / 339200.0f * k4 +
		187.0f / 2100.0f 	* k5 +
		1.0f / 40.0f 		* k6 ;

	v_rk5 = \
		35.0f / 384.0f 		* k0 +
		//0.0f 				* k1 +
		500.0f / 1113.0f 	* k2 +
		125.0f / 192.0f 	* k3 -
		2187.0f / 6784.0f 	* k4 +
		11.0f / 84.0f 		* k5 ;

};

void TimeStepControllerCPU::CalculateDopri5Integration_ATS
(
	Particle& p, Vector& v,
	float timestep, float a, float ah,
	int ts_a, int ts_b
)
{

	float t = 0.0f; // current time during the adaptive timestep estimation
	float T = timestep;  // maximum time for the current advection step

	// the dynamically adjusted current stepsize
	float h = timestep;

	// the conversion factor between relative and real timestep
	float f = ah / h;

	// set current stepsize to last one, if inside bounds
	if(( p.Stepsize>0 ) && ( p.Stepsize<h )) {
		h = p.Stepsize;
	} else {
		p.Stepsize = h;
	}

	Vector pos = p.Position;

	while( !p.bOutOfField && h>0 )
	{
		Vector velocity_rk4, velocity_rk5;

		// Calculate 4-th and 5-th order integration with Dopri-5 integration scheme
		CalculateDopri5Integration(p, velocity_rk4, velocity_rk5, h, a, f*h, ts_a, ts_b);

		// calculate of error
		float error = VectorLength( velocity_rk5 - velocity_rk4 );

		// Accept current stepsize if error is
		// below threshold or minimum stepsize reached.
		if( ((error / h ) <= m_ats_params.tol ) || ( h <= m_ats_params.h_min ))
		{
			// set new position to evaluate
			pos += h * velocity_rk5;

			// increment time and relative time
			t += h;
			a += (f*h);

			// estimate new stepsize
			h = MAX( m_ats_params.h_min, MIN( m_ats_params.eta * h, m_ats_params.rho * pow((m_ats_params.tol/error * pow(h,5)),0.25f)));

			if(( t+h ) > T ) {
				h = T-t; // The remaining timestep, typically much smaller than h
			} else {
				p.Stepsize = h;
			}
		}
		// Otherwise the stepsize is halved
		else
		{
			h = h/2.0f;
		}
	}

	// set velocity to the difference between previous and the next position, divided by the stepsize
	v = (pos - p.Position) / timestep;
};
