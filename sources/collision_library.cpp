#include "collision_library.h"

#include <chrono>
using namespace std::chrono_literals;


namespace collision
{


//sphere vs cylinder collision detection
CollisionState detectCollision ( DynamicPSphere&  S,
                                 StaticPCylinder& C, seconds_type dt){
    auto &_cy = const_cast<StaticPCylinder &> (C);
    const auto tMin= S.curr_t_in_dt;
    const auto tMax = dt;
    const auto newDT = tMax - tMin;
    const auto S_pos = S.getMatrixToScene() * S.getPos();
    const auto _sradius = S.getRadius();
    const auto cy_pos = _cy.evaluateParent(0.5f,0.5f,1,1);
    const auto _cyradius = _cy.getRadiusX();
    const auto ds = S.computeTrajectory(newDT);
    const auto u = cy_pos(1)(0);
    const auto v = cy_pos(0)(1);
    const auto n = u ^ v;
    const auto _n = GMlib::Vector<float,3>(n).getNormalized();
    const auto epsilon = 0.000001;
    const auto _d = (cy_pos(0)(0) + (_sradius + _cyradius) *_n) - S_pos;
    const auto _R = ((S.computeTrajectory(newDT)* _n));
    const auto _Q =(_d * _n);
    if (std::abs(_Q) < epsilon){
        return CollisionState(seconds_type(0.0),CollisionStateFlag::SingularityParallelAndTouching);
    }
    else if(std::abs(_R) < epsilon){
        return CollisionState(seconds_type(0.0),CollisionStateFlag::SingularityParallel);
    }

    const auto x = _Q / _R;
    return CollisionState((x*newDT)+tMin);
}


void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>&        S,
                        StaticPhysObject<GMlib::PCylinder<float>>& C,
                       seconds_type                                     dt){
    auto &_cy=const_cast<StaticPhysObject<GMlib::PCylinder<float>>&>(C);
    const auto cylinder_pos= _cy.evaluateParent(0.5f,0.5f,1,1);
    const auto u = cylinder_pos(1)(0);
    const auto v = cylinder_pos(0)(1);
    const auto n = u ^ v;
    const auto _n = GMlib::Vector<float,3>(n).getNormalized();
    const auto _currentV = S.velocity;
    const auto _vn = _currentV * _n;
    const auto _newV = _currentV -((2*_vn )*_n);
    S.velocity = _newV;
}



CollisionState
detectCollision ( DynamicPhysObject<GMlib::PSphere<float>>& S0,
                  DynamicPhysObject<GMlib::PSphere<float>>& S1,
                 seconds_type                                    dt)
{
    const auto dt_max = dt;
    const auto dt_min = std::max(S0.curr_t_in_dt, S1.curr_t_in_dt);
    const auto new_dt = dt_max - dt_min;

    const auto S0_position = S0.getMatrixToScene() * S0.getPos();
    const auto S1_position = S1.getMatrixToScene() * S1.getPos();

    const auto S0_radius = S0.getRadius();
    const auto S1_radius = S1.getRadius();
    const auto radius_sum = S0_radius + S1_radius;

    const auto Q = (S1_position - S0_position);
    auto r1 = S1.computeTrajectory(new_dt);
    auto r2 = S0.computeTrajectory(new_dt);


    if( S1._state == States::Rolling ) {
        auto unconst_S1 = const_cast<DynamicPhysObject<GMlib::PSphere<float>>&>(S1);
        r1 = unconst_S1.adjustedTrajectory(new_dt);
    }

    if( S0._state == States::Rolling ) {
        auto unconst_S0 = const_cast<DynamicPhysObject<GMlib::PSphere<float>>&>(S0);
        r2 = unconst_S0.adjustedTrajectory(new_dt);
    }

    const auto R = r1 - r2;

    const auto _QR = Q * R;
    const auto _QRQR = std::pow( _QR, 2);

    const auto _RR = R * R;
    const auto _QQ = Q * Q;

    const auto _rr = std::pow( radius_sum, 2);
    const auto _square = std::sqrt(_QRQR - (_RR * (_QQ - _rr)));

    const auto epsilon = 0.00001;

    if ( _square < 0 )
    {
        return CollisionState(seconds_type(0.0), CollisionStateFlag::SingularityNoCollision);
    }
    else if ( (_QQ - _rr) < epsilon )
    {
        return CollisionState(seconds_type(0.0), CollisionStateFlag::SingularityParallelAndTouching);
    }
    else if ( _RR < epsilon )
    {
        return CollisionState( seconds_type(0.0), CollisionStateFlag::SingularityParallel);
    }

    const auto x = (-_QR - _square) / _RR;

    return CollisionState(((x * new_dt) + dt_min), CollisionStateFlag::Collision);


}

void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>& S0,
                       DynamicPhysObject<GMlib::PSphere<float>>& S1,
                       seconds_type                              dt)
{
    const auto S0_old_vel = S0.velocity;    // 2.1
    const auto S1_old_vel = S1.velocity;    // -2.1

    const auto S0_pos = S0.getPos().toType<double>();
    const auto S1_pos = S1.getPos().toType<double>();

    const auto S0_mass = S0.mass;
    const auto S1_mass = S1.mass;

    const auto distance_vector_d = GMlib::Vector<double,3>(S1_pos - S0_pos);
    const auto normal_d = distance_vector_d.getNormalized();
    const auto n = (GMlib::Vector<double,3>(distance_vector_d).getLinIndVec()).getNormalized();

    const auto v0_d = (S0_old_vel * normal_d);
    const auto v1_d = (S1_old_vel * normal_d);
    const auto v0_n = (S0_old_vel * n);
    const auto v1_n = (S1_old_vel * n);

    const auto new_v0_d = (((S0_mass - S1_mass) / (S0_mass + S1_mass) ) * v0_d ) + (((2 * S1_mass) / (S0_mass + S1_mass) ) * v1_d );
    const auto new_v1_d = (((S1_mass - S0_mass) / (S0_mass + S1_mass) ) * v1_d ) + (((2 * S0_mass) / (S0_mass + S1_mass) ) * v0_d );

    const auto S0_new_vel = (v0_n * n) + (new_v0_d * normal_d);  // -2.1
    const auto S1_new_vel = v1_n * n + new_v1_d * normal_d;     // 2.1

    S0.velocity = S0_new_vel;
    S1.velocity = S1_new_vel;
}

CollisionState
detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S,
                 const StaticPhysObject<GMlib::PPlane<float>>&   P,
                 seconds_type                                    dt)
{
    const auto dt_max = dt;
    const auto dt_min = S.curr_t_in_dt;
    const auto new_dt = dt_max - dt_min;

    const auto s_position = S.getMatrixToScene() * S.getPos();
    const auto s_radius = S.getRadius();

    auto &plane = const_cast<StaticPhysObject<GMlib::PPlane<float>>&>(P);
    const auto p = plane.evaluateParent(0.5f, 0.5f, 1, 1);                      // plane.getMatrixToScene * plane.evaluateParent(0.5f, 0.5f, 1, 1)
    const auto plane_pos = p(0)(0);
    const auto u = p(1)(0);
    const auto v = p(0)(1);


    const auto n = u ^ v;
    const auto n_normal = GMlib::Vector<float,3>(n).getNormalized();

    const auto d = (plane_pos + s_radius * n_normal) - s_position;

    auto ds = S.computeTrajectory(new_dt);

    // If the sphere's state is Rolling, it should use an adjusted DS
    if( S._state == States::Rolling ) {
        auto unconst_S = const_cast<DynamicPhysObject<GMlib::PSphere<float>>&>(S);
        ds = unconst_S.adjustedTrajectory(new_dt);
    }

    const auto Q = (d * n_normal);
    const auto R = ( ds * n_normal );    // S.computeTrajectory(dt)

    const auto epsilon = 0.00001;

    if( std::abs(Q) < epsilon )
    {
        // The sphere is "touching" the surface
        return CollisionState( seconds_type(0.0), CollisionStateFlag::SingularityParallelAndTouching);
    }
    else if( std::abs(R) < epsilon)
    {
        return CollisionState( seconds_type(0.0), CollisionStateFlag::SingularityParallel);
    }

    const auto x = Q / R;

    return CollisionState( (x * new_dt) + dt_min);
}

void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>&      S,
                       const StaticPhysObject<GMlib::PPlane<float>>&  P,
                       seconds_type                                   dt) {


    auto &unconst_P = const_cast<StaticPhysObject<GMlib::PPlane<float>>&>(P);
    auto s_pos =S.getPos();
    auto _Radius = S.getRadius();
    auto plane_pos = unconst_P.evaluateParent(0.5f, 0.5f,1, 1);

    auto v = plane_pos(0)(1);
    auto u = plane_pos(1)(0);
    auto n = u^v;
    auto nNormal = GMlib::Vector<float,3>(n).getNormalized();

    auto newVel = (S.velocity - 2*(S.velocity*nNormal)*nNormal)*0.95;
    S.velocity = newVel;


}

CollisionState detectCollision (const DynamicPSphere&  S,
                                const StaticPBezierSurf& B, seconds_type dt) {

    const auto dt_max = dt;
    const auto dt_min = S.curr_t_in_dt;
    const auto new_dt = dt_max - dt_min;

    const auto p0 = S.getPos();
    const auto r = S.getRadius();

    float u, v, t;
    u = 0.5;
    v = 0.5;
    t = 0.0;
    const auto epsilon = 1e-5;

    for ( int i = 0; i < 50; i++) {

        auto ds = S.computeTrajectory( new_dt );
        auto &Surf = const_cast<StaticPBezierSurf&>(B);
        const auto surf = Surf.evaluate(u, v, 1, 1);
        const auto p = p0 + ds*t;
        const auto q = surf(0)(0);
        const auto Su = surf(1)(0);
        const auto Sv = surf(0)(1);
        const auto Sn = GMlib::Vector<float,3>(Su ^ Sv).getNormalized();



        GMlib::SqMatrix<float,3> A;
        A.setCol(Su, 0);
        A.setCol(Sv, 1);
        A.setCol(-ds, 2);
        auto A_inv = A;

        A_inv.invert();

        const auto b = GMlib::Vector<float,3> {p-q-Sn*r};

        const auto x = A_inv * b;

        const auto deltaU = x(0);
        const auto deltaV = x(1);
        const auto deltaT = x(2);

        u += deltaU;
        v += deltaV;
        t += deltaT;

        if( (std::abs(deltaU) < epsilon) and (std::abs(deltaV) < epsilon) and (std::abs(deltaT) < epsilon) ) {

            return CollisionState(seconds_type(deltaT), CollisionStateFlag::Collision);

        }
    }

    return CollisionState(seconds_type(dt_min), CollisionStateFlag::SingularityNoCollision);


}
//Bjorn's code

GMlib::Vector<float,3>
MyController::closestPoint(DynamicPSphere* sphere, seconds_type dt)
{
//    float u, v;
//    P.estimateClpPar(S.getPos(), u, v);
//    auto closest_p = P.getClosestPoint(S.getPos(),u,v);
//    return closest_p;
    auto max_dt = dt;
    auto min_dt = sphere->curr_t_in_dt;
    auto new_dt = max_dt -min_dt;
    auto epsilon= 1e-5;
    float u = 0.4;
    float v = 0.4;
    float delta_u = 0.5;
    float delta_v = 0.5;
    auto  s = sphere->getMatrixToScene() * sphere->getPos();
    auto r = sphere->getRadius();
    GMlib::Vector<double, 3> ds = sphere->computeTrajectory(new_dt);
    auto p = s+ds;
    GMlib::SqMatrix<float,2> A;
    GMlib::Vector<float, 2> b;
    GMlib::Vector<float,3> d{0.0f,0.0f,0.0f};


    auto  planes = _attachedPlanes[sphere];


    //iteration

    for ( int i=0; i<10;i++/* delta_u > epsilon && delta_v > epsilon*/){
    GMlib::Vector <float,3>Sn {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Su {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Sv {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Suu {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Svv {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Suv {0.0f,0.0f,0.0f};
    GMlib::APoint <float,3>q;

    for(auto it = planes.begin(); it != planes.end(); it++){
        GMlib::Vector<float,3> normal{0.0f,0.0f,0.0f};
        GMlib::DMatrix<GMlib::Vector<float,3>>  M = (*it)->evaluateParent(u,v,2,2);
         q     = M(0)(0);
         Su += M(1)(0);
         Sv += M(0)(1);
        Suu += M(2)(0);
        Svv += M(0)(2);
        Suv += M(1)(1);
   }
    d =(q-p);
    A[0][0] = d* Suu + Su * Su;
    A[0][1] = d* Suv + Su * Sv;
    A[1][0] = d* Suv + Su * Sv;
    A[1][1] = d* Svv + Sv * Sv;
    GMlib::SqMatrix<float,2> A_inv = A;
    A_inv.invert();
    b[0] = - d * Su;
    b[1] = - d * Sv;
    GMlib::APoint<float, 3> X = A_inv*b;
    delta_u = X(0);
    delta_v = X(1);
    u += delta_u;
    v += delta_v;


}
    return d;
}


CollisionState detectCollision (const DynamicPSphere& S0,
                                const StaticPSphere& S1, seconds_type dt) {

}

void
computeImpactResponse (DynamicPSphere& S0, const StaticPSphere& S1,
                       seconds_type dt) {

}




void
DynamicPhysObject<GMlib::PSphere<float> >::simulateToTInDt(seconds_type t){

    if( this->_state == States::Still
            || (this->_state == States::Rolling and std::abs(this->velocity(2) <= 1e-2))) {

        //this->curr_t_in_dt = t;
        this->velocity = {0.0f, 0.0f, 0.0f};
        this->_state = States::Still;
        this->environment = &this->_sphereController->noGravEnv;
        //this->rotateGlobal(GMlib::Angle(ds.getLength()/this->getRadius()),GMlib::Vector<float,3>{0,0,1}^ds);

    }
    else {

        auto dt = (t - this->curr_t_in_dt);
        auto MI = this->getMatrixToSceneInverse();

        //move
        auto planes = this->_sphereController->getAttachedObjects (this);
        GMlib::Vector <float,3>n {0.0f,0.0f,0.0f};
        for (auto &it :planes){
            auto M = it->evaluateParent(0.5f,0.5f,1,1);
            auto u = M(1)(0);
            auto v = M(0)(1);
            auto normal = GMlib::Vector<float,3>(u ^ v);
            n+=normal;
           }
        n= GMlib::Vector <float,3>(n/planes.size()).getNormalized();

        GMlib::Vector<double,3> ds {0.0f, 0.0f, 0.0f};

        if( this->_state == States::Rolling ) {

            ds = adjustedTrajectory(dt);
            this->environment = &this->_sphereController->_env;
        }
        if( this->_state == States::Free ) {

            ds = this->computeTrajectory(dt);
            this->environment = &this->_sphereController->_env;
        }
        else
            ds = this->computeTrajectory(dt);

        // Move
        this->translateParent( MI * ds);
        this->curr_t_in_dt = t;

        //std::cout << ds.getLength() <<"  "<< this->getRadius() <<"  "<< (GMlib::Vector<float,3>(0,0,1)^ds) << std::endl;

        if ((ds.getLength() > 1e-2)&& this->_state != States::Still)
        this->rotate(GMlib::Angle(ds.getLength()/this->getRadius()),GMlib::Vector<float,3>(0,0,1)^ds);

        // Update physics
        auto F = this->externalForces();
        auto c = dt.count();
        auto a = F*c;
        this->velocity += a;

        //this->environment = &this->_sphereController->_env;

//        std::cout << this->getPos()(0) <<"  "<< this->getPos()(1) <<"  "<< this->getPos()(2) << std::endl;
    }

}


GMlib::Vector<double,3>
DynamicPhysObject<GMlib::PSphere<float> >::computeTrajectory(seconds_type dt) const {

    auto vel = this->velocity;
    auto dtCount = dt.count();  // 0
    //auto ds = vel * dtCount + 0.5 * xF * std::pow(dtCount, 2);
    GMlib::Vector<double,3> ds = vel * dtCount + 0.5 * this->externalForces() * std::pow(dtCount, 2);

    return ds;

}



GMlib::Vector<double,3>
DynamicPhysObject<GMlib::PSphere<float> >::adjustedTrajectory(seconds_type dt) {


    // Update ds to modified DS'
    auto ds = this->computeTrajectory(seconds_type(dt));
    auto r = this->getRadius();
    auto s = this->getMatrixToScene() * this->getPos();
    //auto p = ds + s;
    auto planes = this->_sphereController->getAttachedObjects(this);
    GMlib::Vector<float,3> n {0.0f, 0.0f, 0.0f};

    //GMlib::Point<float,3> q;

    for( auto& plane : planes ) {

        const auto M = plane->evaluateParent(0.5f, 0.5f, 1, 1);
        const auto u = M(1)(0);
        const auto v = M(0)(1);
        auto normal = GMlib::Vector<float,3> (u ^ v);
        n += normal;
        //q = _sphereController->closestPoint(*this, dt);
    }

    n = GMlib::Vector<float,3> (n / planes.size()).getNormalized();
    auto closest_p = this->_sphereController->closestPoint(this,dt);

    auto adjusted_ds = ds+ (n*r)+closest_p;

    //std::cout << "Good" << std::endl;

    return adjusted_ds;

}

GMlib::Vector<double,3>
DynamicPhysObject<GMlib::PSphere<float> >::externalForces() const {

    assert(environment != nullptr);

    return this->environment->externalForces().toType<double>();
}

void
MyController::localSimulate(double dt) {


    // Reset time variable for all objects
    for( auto sphere : _dynamic_spheres) {

        sphere->curr_t_in_dt = seconds_type{0.0};
    }

    // Detect state changes and fill up our state container
    detectStateChanges(dt);
    sortAndMakeUniqueStates(_singularities);

    // Collision detection algorithm
    for( auto& sphere : _dynamic_spheres) {

        // Sending in sphere twice in the initial call because the function will require a second object in future calls
        sphereDynamicCollision(sphere, sphere, seconds_type(dt));

        for( auto& plane : _static_planes ) {

            sphereStaticCollision(sphere, plane, seconds_type(dt));
        }
    }
    // Make Collision unique
    sortAndMakeUnique(_collisions);

    // Make both collisions and states unique in relation to each other
    if( !_collisions.empty() and !_singularities.empty() ) {

        crossUnique(_collisions, _singularities);
    }
    else {
        // Make sure that the newest event is at the front of the vector
        std::reverse(_singularities.begin(), _singularities.end() );
        std::reverse(_collisions.begin(), _collisions.end());
    }

    while( !_collisions.empty() or !_singularities.empty() ) {

        // IF BOTH NOT EMPTY
        if( !_collisions.empty() and !_singularities.empty() ) {

            const auto col_time = _collisions.back().t_in_dt;
            const auto sing_time = _singularities.back().time;

            // Resolve Collision
            if( col_time < sing_time ) {

                auto c = _collisions.back();
                _collisions.pop_back();;

                detectCollisions(c, dt);     // Also detects more collisions

                detectStateChanges(dt);

                sortAndMakeUnique(_collisions);
                sortAndMakeUniqueStates(_singularities);

                if( !_collisions.empty() and !_singularities.empty() ) {

                    crossUnique(_collisions, _singularities);
                }
                else {
                    // Make sure that the newest event is at the front of the vector
                    std::reverse(_singularities.begin(), _singularities.end() );
                    std::reverse(_collisions.begin(), _collisions.end());
                }
            }

            // Resolve Singularity
            else {

                auto s = _singularities.back();
                _singularities.pop_back();

                detectStates(s, dt);

                // Collision detection algorithm
                for( auto& sphere : _dynamic_spheres) {

                    // Sending in sphere twice in the initial call because the function will require a second object in future calls
                    sphereDynamicCollision(sphere, sphere, seconds_type(dt));

                    for( auto& plane : _static_planes ) {

                        sphereStaticCollision(sphere, plane, seconds_type(dt));
                    }
                }

                detectStateChanges(dt);

                sortAndMakeUnique(_collisions);
                sortAndMakeUniqueStates(_singularities);

                if( !_collisions.empty() and !_singularities.empty() ) {

                    crossUnique(_collisions, _singularities);
                }
                else {
                    // Make sure that the newest event is at the front of the vector
                    std::reverse(_singularities.begin(), _singularities.end() );
                    std::reverse(_collisions.begin(), _collisions.end());
                }
            }
        }

        // IF COLLISIONS NOT EMPTY
        else if( !_collisions.empty() and _singularities.empty() ) {

            auto c = _collisions.back();
            _collisions.pop_back();;

            detectCollisions(c, dt);     // Also detects more collisions

            detectStateChanges(dt);

            sortAndMakeUnique(_collisions);
            sortAndMakeUniqueStates(_singularities);

            if( !_collisions.empty() and !_singularities.empty() ) {

                crossUnique(_collisions, _singularities);
            }
            else {
                // Make sure that the newest event is at the front of the vector
                std::reverse(_singularities.begin(), _singularities.end() );
                std::reverse(_collisions.begin(), _collisions.end());
            }

        }

        // IF SINGULARITIES NOT EMPTY
        else if( _collisions.empty() and !_singularities.empty() ) {

            auto s = _singularities.back();
            _singularities.pop_back();

            detectStates(s, dt);

            // Collision detection algorithm
            for( auto& sphere : _dynamic_spheres) {

                // Sending in sphere twice in the initial call because the function will require a second object in future calls
                sphereDynamicCollision(sphere, sphere, seconds_type(dt));

                for( auto& plane : _static_planes ) {

                    sphereStaticCollision(sphere, plane, seconds_type(dt));
                }
            }

            detectStateChanges(dt);

            sortAndMakeUnique(_collisions);
            sortAndMakeUniqueStates(_singularities);

            if( !_collisions.empty() and !_singularities.empty() ) {

                crossUnique(_collisions, _singularities);
            }
            else {
                // Make sure that the newest event is at the front of the vector
                std::reverse(_singularities.begin(), _singularities.end() );
                std::reverse(_collisions.begin(), _collisions.end());
            }
        }
    }


    //         Start simulation for all objects
    for( auto sphere : _dynamic_spheres) {

        sphere->simulateToTInDt(seconds_type(dt));
    }

}


// Singularity handeling
void
MyController::detectStates(StateChangeObj &state, double dt) {

    auto sphere = state.obj;
    auto newState = state.state;
    auto time = state.time;
    auto planes = state.attachedPlanes;

    if( newState == States::Free ) {

        detachObjects(sphere);
        sphere->_state = newState;
    }
    else {

        setAttachedObjects(planes, sphere);
        sphere->_state = newState;

        if( newState == States::Still) {

            sphere->environment = &noGravEnv;
            sphere->velocity = GMlib::Vector<double,3> (0.0f, 0.0f, 0.0f);
        }
    }

    sphere->simulateToTInDt(time);
}

// Collision handeling
void
MyController::detectCollisions(CollisionObject &c, double dt) {


    // Add more objects here if you end up using more
    auto d_sphere_1     = dynamic_cast<DynamicPSphere*>(c.obj1);
    auto d_sphere_2     = dynamic_cast<DynamicPSphere*>(c.obj2);
    auto s_sphere_2     = dynamic_cast<StaticPSphere*>(c.obj2);
    auto s_plane_2      = dynamic_cast<StaticPPlane*>(c.obj2);
   //auto s_bezier_2     = dynamic_cast<StaticPBezierSurf*>(c.obj2);



    // Impact response
    // If the first object is a sphere
    if(d_sphere_1) {


        if (d_sphere_2)
            collision::computeImpactResponse( *d_sphere_1, *d_sphere_2, c.t_in_dt);    // D_Sphere vs. D_Sphere
        else if (d_sphere_1 && s_sphere_2)
            collision::computeImpactResponse( *d_sphere_1, *s_sphere_2, c.t_in_dt);    // D_Sphere vs. S_Sphere

        else if (d_sphere_1 && s_plane_2) {

            if(d_sphere_1->_state != States::Still)
            {
                d_sphere_1->simulateToTInDt(c.t_in_dt);
            }
            collision::computeImpactResponse( *d_sphere_1, *s_plane_2, c.t_in_dt);      // D_Sphere vs. S_Plane

        }

    }

    // If the dynamic object (obj1) is a sphere
    if( d_sphere_1) {

        sphereDynamicCollision(d_sphere_1, c.obj2, seconds_type(dt));   // Does it collide with any dynamic objects? Can't with same obj as last time
        sphereStaticCollision(d_sphere_1, c.obj2, seconds_type(dt));    // Does it collide with any static objects?  Can't with same obj as last time

        // If sphere 1 collided with a dynamic sphere, check for that sphere's future collisions
        if(d_sphere_2) {

            sphereDynamicCollision(d_sphere_2, d_sphere_1, seconds_type(dt));   // Does it collide with any dynamic objects? Can't with sphere 1
            sphereStaticCollision(d_sphere_2, c.obj2, seconds_type(dt));        // Does it collide with any static objects? Placeholder variable for "Last object"
        }


        if(d_sphere_2) {

            sphereDynamicCollision(d_sphere_2, d_sphere_1, seconds_type(dt));
            sphereStaticCollision(d_sphere_2, c.obj2, seconds_type(dt));
        }
    }
}



void
MyController::detectStateChanges(double dt) {

    for( auto& sphere : _dynamic_spheres) {

        auto singularity = detectStateChange(sphere, dt);

        if (singularity.state != sphere->_state) {

            _singularities.push_back(singularity);
        }
    }

}

//**** Code developed with help from Ghada Bouzidi and Fatemeh Heidari *****

StateChangeObj
MyController::detectStateChange(DynamicPSphere *sphere, double dt) {

    std::unordered_set<StaticPPlane*> planeContainer;       // Used for returning planes that the sphere is (not) attached to
    States state;                           // Holder variable for the state the sphere will enter

    const auto epsilon = 1e-5;

    // Sphere variables
    const auto r = sphere->getRadius();
    const auto pos = sphere->getMatrixToScene() * sphere->getPos();

    // Time variables
    const auto sphereTime = sphere->curr_t_in_dt;
    const auto maxDt = seconds_type(dt);
    const auto newDt = maxDt - sphereTime;
    seconds_type returnTime = sphereTime;
    //        auto still = 0;                                   // Needed to print out the still value further down

    const auto ds = sphere->computeTrajectory(newDt);   // Calculating "original ds"


    // Plane variables.
    auto planes = getAttachedObjects(sphere);
    GMlib::APoint<float,3> q;
    GMlib::Vector<float,3> n {0.0f, 0.0f, 0.0f};

    if( planes.empty() ) {  // Sphere NOT attached

        for (auto& plane : _static_planes) {
            auto M = plane->evaluateParent(0.5f,0.5f,1,1);
            auto q = M(0)(0);
            auto u = M(1)(0);
            auto v = M(0)(1);
            auto n = GMlib::Vector<float,3>(u ^ v).getNormalized();
            auto d = (q + r * n) - pos;

            auto bla        = std::abs(((-n*r) * ds) -(ds*ds));
            auto dsn        = ds * n;
            auto dn         = d*n;
             auto x    = dn / dsn;
            if (x == 0 || x == -INFINITY){
                x = sphereTime.count();
            }
            returnTime      = (x * newDt) + sphereTime;


            if( std::abs(dn) < epsilon and dsn <= 0 ) {

                planeContainer.insert(plane);
                state = States::Rolling;
                std::cout << "detectStateChange says the state will become Rolling from Free" << std::endl;
            }
            else if(std::abs(dn) < epsilon and bla < epsilon ) {

                planeContainer.insert(plane);
                state = States::Still;
                std::cout << "detectStateChange says the state will become Still from Free" << std::endl;
            }
            else state = States::Free;
        }

        return StateChangeObj(sphere, planeContainer, returnTime, state);
    }
    else {      // Sphere ATTACHED

        for (auto &it :planes){
            auto M = it->evaluateParent(0.5f,0.5f,1,1);
            auto pos= M(0)(0);
            auto u = M(1)(0);
            auto v = M(0)(1);
            auto normal = GMlib::Vector<float,3>(u ^ v);
            n+=normal;
            q=pos;
        }
        n= GMlib::Vector <float,3>(n/planes.size()).getNormalized();

        auto d       = (q + r * n) - pos;
        auto bla     = std::abs(((-n*r) * ds) -(ds*ds));
        auto dsn     = ds * n;
        auto dn      = d * n;
        auto x = dn / dsn;
        if (x == 0 || x == -INFINITY){
            x = sphereTime.count();
        }
        returnTime   = (x * newDt) + sphereTime;


        if( sphere->_state == States::Rolling ) {

            if( dsn > 0) {

                std::cout << "detectStateChange says the state will become Free from Rolling" << std::endl;
                state = States::Free;
                return StateChangeObj(sphere, planes, returnTime, state);
            }
            else if( bla < epsilon ) {

                std::cout << "detectStateChange says the state will become Still from Rolling" << std::endl;
                state = States::Still;
                return StateChangeObj(sphere, planes, returnTime, state);
            }
            else return StateChangeObj(sphere, planes, returnTime,States::Rolling);
        }

        else if( sphere->_state == States::Still ) {

            if( bla > epsilon ) {

                std::cout << "detectStateChange says the state will become Rolling from Still" << std::endl;
                state = States::Rolling;
                return StateChangeObj(sphere, planes, returnTime, States::Rolling);
            }
            else if( dsn > 0) {
                std::cout << "detectStateChange says the state will become Free from Still" << std::endl;
                state = States::Free;
                return StateChangeObj(sphere, planes, returnTime, States::Rolling);
            }
            else return StateChangeObj(sphere, planes, returnTime, States::Still);
        }
    }

}


// Get objects attached to sphere
std::unordered_set<StaticPPlane *>
MyController::getAttachedObjects(DynamicPSphere* sphere)
{
    static std::unordered_set<StaticPPlane*> empty {};
    auto iter =_attachedPlanes.find(sphere);

    if( iter !=_attachedPlanes.end() ) {

        return iter->second;
    }
    else return empty;
}

// Set objects attached to sphere
void
MyController::setAttachedObjects(std::unordered_set<StaticPPlane *> plane, DynamicPSphere* sphere)
{
    for( auto& p : plane) {
       _attachedPlanes[sphere].emplace(p);
    }
}

// Remove objects from the set In the map
void
MyController::detachObjects(DynamicPSphere *sphere){

   _attachedPlanes.erase(sphere);

}

// Adding objects to vector

void
MyController::add(DynamicPSphere * const sphere) {

    sphere->environment = &_env;
    _dynamic_spheres.push_back(sphere);
   _attachedPlanes[sphere];

}
void
MyController::add(StaticPSphere * const sphere) {

    _static_spheres.push_back(sphere);
}
void
MyController::add(StaticPPlane * const plane) {

    _static_planes.push_back(plane);
}
void
MyController::add(StaticPCylinder * const cylinder) {

    _static_cylinders.push_back(cylinder);
}
void
MyController::add(StaticPBezierSurf * const surf) {

    _static_bezier_surf.push_back(surf);
}


void collision::DynamicPhysObject<GMlib::PSphere<float> >::getThrough()
{
    if (this->going_through == false){
        this->going_through = true;
        //this->environment = &this->_sphereController->noGravEnv;
        std::cout << "true" << std::endl;

    }
    else if (this->going_through == true){
        this->going_through = false;
       // this->environment = &this->_sphereController->_env;
        std::cout << "false" << std::endl;
    }

    //this->translateGlobal(GMlib::Vector<float,3>(2.0f,0.0f,0.0f));

}




void collision::DynamicPhysObject<GMlib::PSphere<float> >::moveUp()
{
    //std::cout << "up" << std::endl;
    this->checker = true;

    if (this->_state == States::Still){
//        this->_state = States::Rolling;

        this->velocity[2] += 1.0;

        //std::cout << "Foo" << std::endl;
    }

    GMlib::Vector<double,3> newVelVect = this->getVelocity();

            if (newVelVect[1] < 8.0 && newVelVect[1] > -8.0)
            {
                if (newVelVect[1] < 0.0)
                {
                    newVelVect[1] = 0.0;
                }

                newVelVect[1] += 1.0;
                newVelVect[0] *= 0.5;

                this->setVelocity(newVelVect);

            }
            else
            {
                while (newVelVect[1] >= 8.0 || newVelVect[1] <= -8.0)
                {
                    newVelVect[1] *= 0.9;
                    this->setVelocity(newVelVect);
                }
            }

}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::moveDown()
{
    //std::cout << "down" << std::endl;

    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }


    GMlib::Vector<double,3> newVelVect = this->getVelocity();
    if (newVelVect[1] < 8.0 && newVelVect[1] > -8.0)
    {
        if (newVelVect[1] > 0.0)
        {
            newVelVect[1] = 0.0;
        }

        newVelVect[1] -= 1.0;
        newVelVect[0] *= 0.5;

        this->setVelocity(newVelVect);
    }
    else
    {
        while (newVelVect[1] >= 8.0 || newVelVect[1] <= -8.0)
        {
            newVelVect[1] *= 0.9;
            this->setVelocity(newVelVect);
        }
    }
}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::moveLeft()
{
    //std::cout << "left" << std::endl;

    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }


    GMlib::Vector<double,3> newVelVect = this->getVelocity();
            if (newVelVect[0] < 8.0 && newVelVect[0] > -8.0)
            {
                if (newVelVect[0] > 0.0)
                {
                    newVelVect[0] = 0.0;
                }

                newVelVect[0] -= 1.0;
                newVelVect[1] *= 0.5;

                this->setVelocity(newVelVect);
            }
            else
            {
                while (newVelVect[0] >= 8.0 || newVelVect[0] <= -8.0)
                {
                    newVelVect[0] *= 0.9;
                    this->setVelocity(newVelVect);
                }
            }

}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::moveRight()
{
    //std::cout << "right" << std::endl;

    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }


    GMlib::Vector<double,3> newVelVect = this->getVelocity();
            if (newVelVect[0] < 8.0 && newVelVect[0] > -8.0)
            {
                if (newVelVect[0] < 0.0)
                {
                    newVelVect[0] = 0.0;
                }

                newVelVect[0] += 1.0;
                newVelVect[1] *= 0.5;

                this->setVelocity(newVelVect);
            }
            else
            {
                while (newVelVect[0] >= 8.0 || newVelVect[0] <= -8.0)
                {
                    newVelVect[0] *= 0.9;
                    this->setVelocity(newVelVect);
                }
            }
}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::setVelocity(const GMlib::Vector<double, 3> vel)
{
    this->velocity = vel;
}

GMlib::Vector<double, 3> collision::DynamicPhysObject<GMlib::PSphere<float> >::getVelocity()
{
    return this->velocity;
}

double collision::DynamicPhysObject<GMlib::PSphere<float> >::getMass()
{
    return this->mass;
}

GMlib::Vector<float, 3> collision::DynamicPhysObject<GMlib::PSphere<float> >::getSurfNormal()
{

//auto _surface

}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::computeStep(double dt)
{

}

int collision::StaticPhysObject<GMlib::PPlane<float> >::getId() const
{
    return this->id;
}

void collision::StaticPhysObject<GMlib::PPlane<float> >::setId(int value)
{
    this->id = value;
}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::translateUp()
{
    this->translateGlobal(GMlib::Vector<float,3>(0.0f,2.0f,0.0f));
}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::translateDown()
{
    this->translateGlobal(GMlib::Vector<float,3>(0.0f,-2.0f,0.0f));
}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::translateLeft()
{
    this->translateGlobal(GMlib::Vector<float,3>(-2.0f,0.0f,0.0f));
}

void collision::DynamicPhysObject<GMlib::PSphere<float> >::translateRight()
{
    this->translateGlobal(GMlib::Vector<float,3>(2.0f,0.0f,0.0f));
}







} // END namespace collision


