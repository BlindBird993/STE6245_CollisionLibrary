#include "collision_library.h"

#include <chrono>
using namespace std::chrono_literals;


namespace collision
{

#define detectCollisionSphereSphereMethod {
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
        r1 = S1.adjustedTrajectory(new_dt);
    }

    if( S0._state == States::Rolling ) {
        r2 = S0.adjustedTrajectory(new_dt);
    }

    const auto R = r1 - r2;

    const auto _QR = Q * R;
    const auto _QRQR = std::pow( _QR, 2);

    const auto _RR = R * R;
    const auto _QQ = Q * Q;

    const auto _rr = std::pow( radius_sum, 2);
    const auto _square = std::sqrt(_QRQR - (_RR * (_QQ - _rr)));

    const auto epsilon = 1e-7;

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

#define end1 }

#define ImpactResponseSphereSphereMethod {
void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>& S0,
                       DynamicPhysObject<GMlib::PSphere<float>>& S1,
                       seconds_type                              dt)
{
    const auto S0_old_vel = S0.velocity;
    const auto S1_old_vel = S1.velocity;

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

    const auto S0_new_vel = (v0_n * n) + (new_v0_d * normal_d);
    const auto S1_new_vel = v1_n * n + new_v1_d * normal_d;

    S0.velocity = S0_new_vel;
    S1.velocity = S1_new_vel;
}
#define end1 }


#define detectCollisionSpherePlaneMethod {
CollisionState
detectCollision (DynamicPhysObject<GMlib::PSphere<float>>& S,
                 StaticPhysObject<GMlib::PPlane<float>>&   P,
                 seconds_type                                    dt)
{
    const auto dt_max = dt;
    const auto dt_min = S.curr_t_in_dt;
    const auto new_dt = dt_max - dt_min;

    const auto s_position = S.getMatrixToScene() * S.getPos();
    const auto s_radius = S.getRadius();

    auto &plane = const_cast<StaticPhysObject<GMlib::PPlane<float>>&>(P);
    const auto p = plane.evaluateParent(0.5f, 0.5f, 1, 1);
    const auto plane_pos = p(0)(0);
    const auto u = p(1)(0);
    const auto v = p(0)(1);


    const auto n = u ^ v;
    const auto n_normal = GMlib::Vector<float,3>(n).getNormalized();

    const auto d = (plane_pos + s_radius * n_normal) - s_position;

    auto ds = S.computeTrajectory(new_dt);

    if( S._state == States::Rolling ) {
        ds = S.adjustedTrajectory(new_dt);
    }

    const auto Q = (d * n_normal);
    const auto R = ( ds * n_normal );

    const auto epsilon = 1e-7;

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
#define end1 }

#define ImpactResponseSpherePlaneMethod {
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
#define end1 }


#define detectCollisionSphereBezierMethod {
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
#define end1 }



#define closestPointMethod {
/** Code developed with help from Bjørn-Richard Pedersen and Fatemeh Heidari **/
// check the planes against a sphere and find a closest one
GMlib::Vector<float,3>
MyController::closestPoint(DynamicPSphere* sphere, seconds_type dt)
{
    auto max_dt = dt;
    auto min_dt = sphere->curr_t_in_dt;
    auto new_dt = max_dt -min_dt;
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

    //use taylor expansion
    //iteration
    for ( int i=0; i<10;i++){
    GMlib::Vector <float,3>Su {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Sv {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Suu {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Svv {0.0f,0.0f,0.0f};
    GMlib::Vector <float,3>Suv {0.0f,0.0f,0.0f};
    GMlib::APoint <float,3>q;

    for(auto it = planes.begin(); it != planes.end(); it++){
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
#define end1 }

#define simulateToTInDtMethod {
/** Code developed with help from Bjørn-Richard Pedersen and Fatemeh Heidari **/
void
DynamicPhysObject<GMlib::PSphere<float> >::simulateToTInDt(seconds_type t){

    if( this->_state == States::Still or
             (this->_state == States::Rolling and std::abs(this->velocity(2) <= 1e-2))) {

        this->velocity = {0.0f, 0.0f, 0.0f};
        this->environment = &this->_sphereController->noGravEnv;

    }
    else {

        auto dt = (t - this->curr_t_in_dt);
        auto MI = this->getMatrixToSceneInverse();

        GMlib::Vector<double,3> ds {0.0f, 0.0f, 0.0f};

        if( this->_state == States::Rolling ) {

            ds = adjustedTrajectory(dt);

            this->translateParent(MI*ds);
            this->curr_t_in_dt = t;

            //Update physics for rolling state
            auto F = this->externalForces();
            auto c = dt.count();
            auto a = F * c;
            this->velocity -= a;

            this->environment = &this->_sphereController->env;
        }

        else{

        ds = this->computeTrajectory(dt);

        // Move
        this->translateParent( MI * ds);
        this->curr_t_in_dt = t;

        //basic rotation, vector is invalid
        if ((ds.getLength() > 1e-2)&& this->_state != States::Still)
        this->rotate(GMlib::Angle(ds.getLength()/this->getRadius()),GMlib::Vector<float,3>(0,0,1)^ds);

        // Update physics
        auto F = this->externalForces();
        auto c = dt.count();
        auto a = F*c;
        this->velocity += a;
        }
    }

}
#define end1 }


#define computeTrajectoryMethod {
GMlib::Vector<double,3>
DynamicPhysObject<GMlib::PSphere<float> >::computeTrajectory(seconds_type dt) const {

    auto vel = this->velocity;
    auto dtc = dt.count();
    GMlib::Vector<double,3> ds = vel * dtc + 0.5 * this->externalForces() * std::pow(dtc, 2);

    return ds;

}
#define end1 }

#define adjustedTrajectoryMethod {
/** Code developed with help from Bjørn-Richard Pedersen and Fatemeh Heidari **/
GMlib::Vector<double,3>
DynamicPhysObject<GMlib::PSphere<float> >::adjustedTrajectory(seconds_type dt) {


    // Update ds
    auto ds = this->computeTrajectory(seconds_type(dt));
    auto r = this->getRadius();
    auto s = this->getMatrixToScene() * this->getPos();
    auto planes = this->_sphereController->getAttachedObjects(this);
    GMlib::Vector<float,3> n {0.0f, 0.0f, 0.0f};

    for( auto& plane : planes ) {

        const auto M = plane->evaluateParent(0.5f, 0.5f, 1, 1);
        const auto u = M(1)(0);
        const auto v = M(0)(1);
        auto normal = GMlib::Vector<float,3> (u ^ v);
        n += normal;
    }

    n = GMlib::Vector<float,3> (n / planes.size()).getNormalized();
    auto closest_p = this->_sphereController->closestPoint(this,dt);

    auto adjusted_ds = ds+ (n*r)+closest_p;

    return adjusted_ds;

}
#define end1 }


#define externalForcesMethod {
GMlib::Vector<double,3>
DynamicPhysObject<GMlib::PSphere<float> >::externalForces() const {

    assert(environment != nullptr);

    return this->environment->externalForces().toType<double>();
}
#define end1 }


#define localSimulateMethod {
/** Code developed with help from Bjørn-Richard Pedersen and Fatemeh Heidari **/
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

            sphereCollision(sphere,seconds_type(dt));
    }
    // Make Collision unique
    sortAndMakeUnique(_collisions);

    // Make both collisions and states unique in relation to each other
    if( !_collisions.empty() and !_singularities.empty() ) {

        crossUnique(_collisions, _singularities);
    }
    else {
        // Make sure that the newest event is at the front of the vector
        reverseMethod(_singularities,_collisions);
    }

    while( !_collisions.empty() or !_singularities.empty() ) {

        // If both containers not empty
        if( !_collisions.empty() and !_singularities.empty() ) {

             auto col_time = _collisions.back().t_in_dt;
             auto sing_time = _singularities.back().time;

            // Resolve Collision
            if( col_time < sing_time ) {

                auto c = _collisions.back();//take a first collision object
                _collisions.pop_back();//remove it from vector

                detectCollisions(c, dt);// Also detects more collisions, handle this collision

                detectStateChanges(dt); //detect States

                sortAndMakeUnique(_collisions);
                sortAndMakeUniqueStates(_singularities);

                if( !_collisions.empty() and !_singularities.empty() ) {

                    crossUnique(_collisions, _singularities);// crosscheck
                }
                else {
                    // Make sure that the newest event is at the front of the vector
                     reverseMethod(_singularities,_collisions);;
                }
            }

            // Resolve Singularity
            else {

                auto s = _singularities.back();
                _singularities.pop_back();

                detectStates(s, dt);

                // Collision detection algorithm
                for( auto& sphere : _dynamic_spheres) {


                        sphereCollision(sphere,seconds_type(dt));

                }

                detectStateChanges(dt);

                sortAndMakeUnique(_collisions);
                sortAndMakeUniqueStates(_singularities);

                if( !_collisions.empty() and !_singularities.empty() ) {

                    crossUnique(_collisions, _singularities);
                }
                else {
                    // Make sure that the newest event is at the front of the vector
                     reverseMethod(_singularities,_collisions);
                }
            }
        }

        // If collisions container not empty
        else if( !_collisions.empty() and _singularities.empty() ) {

            auto c = _collisions.back();
            _collisions.pop_back();;

            detectCollisions(c, dt);     // Also detects more collisions

            detectStateChanges(dt);     //detect state changes

            sortAndMakeUnique(_collisions);
            sortAndMakeUniqueStates(_singularities);

            if( !_collisions.empty() and !_singularities.empty() ) {

                crossUnique(_collisions, _singularities);
            }
            else {
                // Make sure that the newest event is at the front of the vector
                reverseMethod(_singularities,_collisions);
            }

        }

        // If singularities container not empty
        else if( _collisions.empty() and !_singularities.empty() ) {

            auto s = _singularities.back();
            _singularities.pop_back();

            detectStates(s, dt); //detect state changes

            // Collision detection algorithm
            for( auto& sphere : _dynamic_spheres) {


                    sphereCollision(sphere, seconds_type(dt));

            }

            detectStateChanges(dt);

            sortAndMakeUnique(_collisions);
            sortAndMakeUniqueStates(_singularities);

            if( !_collisions.empty() and !_singularities.empty() ) {

                crossUnique(_collisions, _singularities);
            }
            else {
                // Make sure that the newest event is at the front of the vector
                reverseMethod(_singularities,_collisions);
            }
        }
    }


    //Start simulation for all objects
    for( auto sphere : _dynamic_spheres) {

        sphere->simulateToTInDt(seconds_type(dt));
    }

}
#define endl1 }


#define detectStatesMethod {
/** Code developed with help from Bjørn-Richard Pedersen **/
// Singularity handeling
void
MyController::detectStates(StateChangeObj &state, double dt) {

    auto sphere = state.obj;
    auto newState = state.state;
    auto time = state.time;
    auto planes = state.planes;

    if( newState == States::Free ) {

        detachObjects(sphere);
        sphere->_state = newState;
    }
    else {

        setAttachedObjects(planes, sphere);
        sphere->_state = newState;
    }

    sphere->simulateToTInDt(time);
}
#define end1 }

#define detectCollisionsMethod {
/** Code developed with help from Bjørn-Richard Pedersen **/
// Collision handeling
void
MyController::detectCollisions(CollisionObject &c, double dt) {// check for collisions inside the template, handle it here


    auto d_sphere_1     = dynamic_cast<DynamicPSphere*>(c.obj1);
    auto d_sphere_2     = dynamic_cast<DynamicPSphere*>(c.obj2);
    auto s_plane_2      = dynamic_cast<StaticPPlane*>(c.obj2);


    // Impact response
    // If the first object is a sphere
    if(d_sphere_1) {

        if (d_sphere_2){

            if (d_sphere_2->_state == States::Still){
                d_sphere_1->simulateToTInDt(c.t_in_dt);

                d_sphere_2->curr_t_in_dt = d_sphere_1->curr_t_in_dt;
                d_sphere_2->environment = &env;
                d_sphere_2->_state = States::Rolling;

                collision::computeImpactResponse( *d_sphere_1, *d_sphere_2, c.t_in_dt);    // D_Sphere vs. D_Sphere
            }

        else {
                d_sphere_1->simulateToTInDt(c.t_in_dt);
                d_sphere_2->simulateToTInDt(c.t_in_dt);

                collision::computeImpactResponse( *d_sphere_1, *d_sphere_2, c.t_in_dt);    // D_Sphere vs. D_Sphere
            }

            }

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

        sphereCollision(d_sphere_1, seconds_type(dt));    // Does it collide with any static objects?  Can't with same obj as last time

        // If sphere 1 collided with a dynamic sphere, check for that sphere's future collisions
        if(d_sphere_2) {

            sphereCollision(d_sphere_2,seconds_type(dt));        // Does it collide with any static objects? Placeholder variable for "Last object"
        }// additional collisions
    }
}
#define end1 }

#define detectStateChangesMethod {
void
MyController::detectStateChanges(double dt) {

    for( auto& sphere : _dynamic_spheres) {

        auto singularity = detectStateChange(sphere, dt);

        if (singularity.state != sphere->_state) {

            _singularities.push_back(singularity);
        }
    }

}
#define end1 }

#define detectStateChangeMethod {
/** Code developed with help from Ghada Bouzidi and Fatemeh Heidari **/
StateChangeObj
MyController::detectStateChange(DynamicPSphere *sphere, double dt) {

    std::unordered_set<StaticPPlane*> planeContainer;       // Used for returning planes that the sphere is (not) attached to
    States state;                           // Holder variable for the state the sphere will enter

    const auto epsilon = 1e-5;

    // Sphere variables
    const auto r = sphere->getRadius();
    const auto pos = sphere->getMatrixToScene() * sphere->getPos();

    // Time variables
    const auto max_dt = seconds_type(dt);
    const auto new_dt = max_dt - sphere->curr_t_in_dt;
    seconds_type min_dt = sphere->curr_t_in_dt;

    const auto ds = sphere->computeTrajectory(new_dt);


    // Plane variables.
    auto planes = getAttachedObjects(sphere);
    GMlib::APoint<float,3> q;
    GMlib::Vector<float,3> n {0.0f, 0.0f, 0.0f};

    if( planes.empty() ) {  // sphere is not attached to plane

        for (auto& plane : _static_planes) {
            auto M = plane->evaluateParent(0.5f,0.5f,1,1);
            auto q = M(0)(0);
            auto u = M(1)(0);
            auto v = M(0)(1);
            auto n = GMlib::Vector<float,3>(u ^ v).getNormalized();

            auto d = (q + r * n) - pos;
            auto dn         = d*n;
            auto x    = dn / (ds * n);

            min_dt      = (x * new_dt) + sphere->curr_t_in_dt;


            if( std::abs(dn) < epsilon and (ds * n) <= 0 ) {

                planeContainer.insert(plane);
                state = States::Rolling;
            }
            else if(std::abs(dn) < epsilon and std::abs(((-n*r) * ds) -(ds*ds)) < epsilon ) {

                planeContainer.insert(plane);
                state = States::Still;
            }
            else state = States::Free;
        }

        return StateChangeObj(sphere, planeContainer, min_dt, state);
    }
    else {      // sphere is attached to plane

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
        auto dn      = d * n;
        auto x = dn / (ds * n);

        min_dt   = (x * new_dt) + sphere->curr_t_in_dt;


        if( sphere->_state == States::Rolling ) {

            if( (ds * n) > 0) {

                state = States::Free;
                return StateChangeObj(sphere, planes, min_dt, state);
            }
            else if( std::abs(((-n*r) * ds) -(ds*ds)) < epsilon ) {

                state = States::Still;
                return StateChangeObj(sphere, planes, min_dt, state);
            }
            else return StateChangeObj(sphere, planes, min_dt,States::Rolling);
        }

        else if( sphere->_state == States::Still ) {

            if( std::abs(((-n*r) * ds) -(ds*ds)) > epsilon ) {

                state = States::Rolling;
                return StateChangeObj(sphere, planes, min_dt, States::Rolling);
            }
            else if( (ds * n) > 0) {
                state = States::Free;
                return StateChangeObj(sphere, planes, min_dt, States::Rolling);
            }
            else return StateChangeObj(sphere, planes, min_dt, States::Still);
        }
    }

}
#define end1 }

#define getAttachObjectsMethod {
// Get planes attached to sphere
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
#define end1 }

#define getAttachObjectsMethod {
// Set objects attached to sphere
void
MyController::setAttachedObjects(std::unordered_set<StaticPPlane *> planes, DynamicPSphere* sphere)
{
    for( auto& plane : planes) {
       _attachedPlanes[sphere].emplace(plane);
    }
}
#define end1 }

#define detachObjectsMethod {
void
MyController::detachObjects(DynamicPSphere *sphere){

   _attachedPlanes.erase(sphere);

}
#define end1 }

#define ContainerReverseMethod {
void MyController::reverseMethod(std::vector<StateChangeObj> singularitiesContainer,std::vector<CollisionObject> collisionsContainer)
{
    std::reverse(singularitiesContainer.begin(), singularitiesContainer.end() );
    std::reverse(collisionsContainer.begin(), collisionsContainer.end());
}
#define end1 }

#define addingMethods {
void
MyController::add(DynamicPSphere * const sphere) {

    sphere->environment = &env;
    _dynamic_spheres.push_back(sphere);

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
#define end1 }


#define movingSphereMethods {

void collision::DynamicPhysObject<GMlib::PSphere<float> >::moveUp()
{


    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }

    GMlib::Vector<double,3> newVelVect = this->getVelocity();

            if (newVelVect[1] < 8.0 and newVelVect[1] > -8.0)
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

    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }


    GMlib::Vector<double,3> newVelVect = this->getVelocity();
    if (newVelVect[1] < 8.0 and newVelVect[1] > -8.0)
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

    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }


    GMlib::Vector<double,3> newVelVect = this->getVelocity();
            if (newVelVect[0] < 8.0 and newVelVect[0] > -8.0)
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

    if (this->_state == States::Still){

        this->velocity[2] += 1.0;

    }


    GMlib::Vector<double,3> newVelVect = this->getVelocity();
            if (newVelVect[0] < 8.0 and newVelVect[0] > -8.0)
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
#define end1 }



#define getAndSetVelocityMethods {
void collision::DynamicPhysObject<GMlib::PSphere<float> >::setVelocity(const GMlib::Vector<double, 3> vel)
{
    this->velocity = vel;
}

GMlib::Vector<double, 3> collision::DynamicPhysObject<GMlib::PSphere<float> >::getVelocity()
{
    return this->velocity;
}
#define end1 }

#define getAndSetIdMethods {
int collision::StaticPhysObject<GMlib::PPlane<float> >::getId() const
{
    return this->id;
}

void collision::StaticPhysObject<GMlib::PPlane<float> >::setId(int value)
{
    this->id = value;
}
#define end1 }

#define changePositionMethods {
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
#define end1 }






} // END namespace collision


