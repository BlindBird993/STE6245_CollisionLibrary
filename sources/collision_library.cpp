#include "collision_library.h"



namespace collision
{
/**

//Bjorn code, had some problems with the old one, need to figure out the case
CollisionState
detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S0,
                 const DynamicPhysObject<GMlib::PSphere<float>>& S1,
                 seconds_type dt)
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
    const auto R = (S1.computeTrajectory(new_dt) - S0.computeTrajectory(new_dt));

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
    //
}


//    CollisionState
//    detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S0,
//                     const DynamicPhysObject<GMlib::PSphere<float>>& S1,
//                     seconds_type                                    dt)
//    {
//        auto dt_max = dt;
//        auto dt_min = std::max(S0.curr_t_in_dt,S1.curr_t_in_dt);//for plane only use sphere
//        auto dt_new = dt_max - dt_min;

//        auto S0_pose = S0.getMatrixToScene()*S0.getPos();
//        auto S1_pose = S1.getMatrixToScene()*S1.getPos();
//        auto S0_rad = S0.getRadius();
//        auto S1_rad = S1.getRadius();

//        auto rad_sum = S0_rad + S1_rad;
//        const auto Q = S1_pose - S0_pose;
//        const auto R = S1.computeTrajectory(dt_new) - S0.computeTrajectory(dt_new);

//        const auto r_square = std::pow(rad_sum,2);
//        const auto x = (-(Q*R) - sqrt(std::pow(Q*R,2) - (R*R)*((Q*Q)-r_square)))/(R*R);
//        auto d = Q + x*R;


//        const auto eps = std::pow(10,-5);

//        //std::cout << eps << std::endl;
//        return CollisionState((x*dt_new)+dt_min);
//        //
//        //auto d_norm




//    }



CollisionState
detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S,
                 const StaticPhysObject<GMlib::PPlane<float>>&   P,
                 seconds_type                                    dt)
{
    //auto unconst_S = const_cast<DynamicPhysObject<GMlib::PSphere<float>>&>(S);

    const auto dt_max = dt;
    const auto dt_min = S.curr_t_in_dt;
    const auto new_dt = dt_max - dt_min;

    auto &unconst_P = const_cast<StaticPhysObject<GMlib::PPlane<float>>&>(P);
    auto s_pos =S.getMatrixToScene()*S.getPos();
    auto _Radius = S.getRadius();
    auto plane_pos = unconst_P.evaluateParent(0.5f, 0.5f,1, 1);

    auto v = plane_pos(0)(1);
    auto u = plane_pos(1)(0);
    auto n = u^v;

    auto nNormal = GMlib::Vector<float,3>(n).getNormalized();
    auto d = (plane_pos(0)(0) + _Radius*nNormal) - s_pos;

    const auto ds = (S.computeTrajectory(new_dt)*nNormal);

    auto x = (d*nNormal)/(S.computeTrajectory(new_dt)*nNormal);

    const auto epsilon = 0.00001;

    if ( std::abs(d*n) < epsilon )
    {
        return CollisionState(seconds_type(0.0), CollisionStateFlag::SingularityParallelAndTouching);
    }
    else if (std::abs(ds)< epsilon )
    {
        return CollisionState( seconds_type(0.0), CollisionStateFlag::SingularityParallel);
    }

    return CollisionState(((x*new_dt)+ dt_min),CollisionStateFlag::Collision);

}

CollisionState
detectCollision (const DynamicPSphere&  S,
                 const StaticPBezierSurf& B, seconds_type dt)
{
    const auto dt_max = dt;
    const auto dt_min = S.curr_t_in_dt;
    const auto new_dt = dt_max - dt_min;

    const auto _Radius = S.getRadius();//r
    const auto s_pos = S.getPos();//p

    float u,v,t;
    t = 0.0;
    u = 0.5;
    v = 0.5;
    const auto epsilon = 1e-5;

    for(int i=0; i<50; i++){

    auto ds = S.computeTrajectory(new_dt);
    auto &unconst_B = const_cast<StaticPBezierSurf&>(B);
    const auto surf_pos = unconst_B.evaluate(u, v,1, 1);

    const auto s_pos_new = s_pos + ds*t;
    const auto q = surf_pos(0)(0);
    const auto Sv = surf_pos(0)(1);
    const auto Su = surf_pos(1)(0);


    //const auto ds_new = ds*t;

    const auto Sn = GMlib::Vector<float,3>(Su ^ Sv).getNormalized();


    GMlib::SqMatrix<float,3>A;
    A.setCol(Su,0);
    A.setCol(Sv,1);
    A.setCol(-ds,2);
    auto A_inv = A;

    A_inv.invert();

    const auto b = GMlib::Vector<float,3>{s_pos_new-q-Sn*_Radius};
    const auto x = A_inv*b;

    const auto deltaU = x(0);
    const auto deltaV = x(1);
    const auto deltaT = x(2);

    u+=deltaU;
    v+=deltaV;
    t+=deltaT;

    if( (std::abs(deltaU) < epsilon) and (std::abs(deltaV) < epsilon) and (std::abs(deltaT) < epsilon) ) {

        return CollisionState(seconds_type(deltaT), CollisionStateFlag::Collision);
    }
}

return CollisionState(seconds_type(dt_min), CollisionStateFlag::SingularityNoCollision);
}

CollisionState
detectCollision (const DynamicPSphere& S, const StaticPCylinder& C,
                 seconds_type dt)
{
    auto &_cy = const_cast<StaticPCylinder &> (C);
    const auto tMin= S.curr_t_in_dt;
    const auto tMax = dt; const auto newDT = tMax - tMin;
    const auto S_pos = S.getMatrixToScene() * S.getPos();
    const auto _sradius = S.getRadius();
    const auto cy_pos = _cy.evaluateParent(0.5f,0.5f,1,1);
    const auto _cyradius = _cy.getRadiusX();

    const auto ds = S.computeTrajectory(newDT);

    const auto u = cy_pos(1)(0);
    const auto v = cy_pos(0)(1);
    const auto n = u ^ v;

    const auto _n = GMlib::Vector<float,3>(n).getNormalized();
    const auto epsilon = 0.000001; const auto _d = (cy_pos(0)(0) + (_sradius + _cyradius) *_n) - S_pos;
    const auto _R = ((S.computeTrajectory(newDT)* _n));
    const auto _Q =(_d * _n); if (std::abs(_Q) < epsilon){

        return CollisionState(seconds_type(0.0),
                              CollisionStateFlag::SingularityParallelAndTouching);
        }
    else if(std::abs(_R) < epsilon){
        return CollisionState(seconds_type(0.0),
                              CollisionStateFlag::SingularityParallel);
            }
    const auto x = _Q / _R;
    return CollisionState((x*newDT)+tMin);
    }

void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>& S0,
                       DynamicPhysObject<GMlib::PSphere<float>>& S1,
                       seconds_type                              dt)

{

    auto s0_pos = S0.getPos().toType<double>();
    auto s1_pos = S1.getPos().toType<double>();
    auto d = (s1_pos - s0_pos);
    auto dNormalized = GMlib::Vector<double,3>(d).getNormalized();

    auto vel0 = S0.velocity;
    auto vel1 = S1.velocity;



    auto n = GMlib::Vector<double,3>(d).getLinIndVec();
    auto nNormal = n.getNormalized();

    auto V0d = (vel0*dNormalized);
    auto V1d = (vel1*dNormalized);

    auto V0n = (vel0*nNormal);
    auto V1n = (vel1*nNormal);

    auto S0mass = S0.mass;
    auto S1mass = S1.mass;

    auto V0d_new = (((S0mass-S1mass)/(S0mass+S1mass))*V0d + ((2*S1mass)/(S0mass+S1mass))*V1d);
    auto V1d_new = (((S1mass-S0mass)/(S0mass+S1mass))*V1d + ((2*S0mass)/(S0mass+S1mass))*V0d);

    S0.velocity = (V0n*nNormal + (V0d_new*dNormalized));
    S1.velocity = (V1n*nNormal + (V1d_new*dNormalized));

}

void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>& S,
                       const StaticPhysObject<GMlib::PCylinder<float>>&   P,
                       seconds_type                              dt)
{
    auto &unconst_P = const_cast<StaticPhysObject<GMlib::PCylinder<float>>&>(P);
    auto s_pos =S.getPos();
    auto _Radius = S.getRadius();
    auto plane_pos = unconst_P.evaluateParent(0.5f, 0.5f,1, 1);

    auto v = plane_pos(0)(1);
    auto u = plane_pos(1)(0);
    auto n = u^v;
    auto nNormal = GMlib::Vector<float,3>(n).getNormalized();

    auto newVel = (S.velocity - 2*(S.velocity*nNormal)*nNormal);
    S.velocity = newVel;

}

void
computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>& S,
                       const StaticPhysObject<GMlib::PPlane<float>>&   P,
                       seconds_type                              dt)
{
    auto &unconst_P = const_cast<StaticPhysObject<GMlib::PPlane<float>>&>(P);
    auto s_pos =S.getPos();
    auto _Radius = S.getRadius();
    auto plane_pos = unconst_P.evaluateParent(0.5f, 0.5f,1, 1);

    auto v = plane_pos(0)(1);
    auto u = plane_pos(1)(0);
    auto n = u^v;
    auto nNormal = GMlib::Vector<float,3>(n).getNormalized();

    auto newVel = (S.velocity - 2*(S.velocity*nNormal)*nNormal);
    S.velocity = newVel;

}

void DynamicPhysObject<GMlib::PSphere<float> >::simulateToTInDt(seconds_type t)
{
    //start
    auto const M = this->getMatrixToSceneInverse();
    auto dt = this->curr_t_in_dt;
    const auto new_dt = t-dt;

    GMlib::Vector<double,3> ds = this->computeTrajectory(new_dt);

    /**
    //correct trajectory, save as a member

//    auto ps = this->getPos() + ds;

//    static auto g = GMlib::Vector<double,3>(0,0,-9.8);
//    StaticPPlane* _surface = this->_plane;

//    _surface->estimateClpPar(this->getPos(),_u,_v);

//    _surface->getClosestPoint(this->getPos() + ds,_u,_v);


//    GMlib::DMatrix<GMlib::Vector<float,3>> sMatrix = _surface->evaluate(_u,_v,1,1);
//    GMlib::UnitVector<float,3> norm = sMatrix[0][1] ^ sMatrix[1][0];

//    ds += sMatrix[0][0] + norm * _radius;
//    this->velocity+=this->curr_t_in_dt.count()*g;
//    this->velocity-=(this->velocity*norm)*norm;

    //move
    this->translateParent(M*ds);
    //physics
    auto F = environment->externalForces();
    auto a = 0.5*F*std::pow(dt.count(),2)*this->mass;
    this->velocity+=a;

}


//GMlib::Vector<double,3> DynamicPhysObject<GMlib::PSphere<float> >::computeTrajectory(seconds_type dt_in) const
//{

//    auto dt=dt_in.count();
//    auto F = environment->externalForces();
//    auto a = 0.5*F*std::pow(dt,2)*this->mass;
//    return this->velocity*dt + a;

//}

GMlib::Vector<double,3> DynamicPhysObject<GMlib::PSphere<float> >::externalForces() const
{

    assert(environment == nullptr);

    return environment->externalForces();

}


void MyController::collisionAlgorithm(seconds_type dt)
{
    std::vector<CollisionObject> C;
    //std::vector<States> S;
    for(auto itr1 = _dynamic_spheres.begin();itr1 != _dynamic_spheres.end()-1;++itr1){
        for(auto itr2 = itr1;itr2 != _dynamic_spheres.end();++itr2){

            auto& sphere1 = *itr1;
            auto& sphere2 = *itr2;
            auto dt_min_Diff = std::max(sphere1->curr_t_in_dt,sphere2->curr_t_in_dt);

            const auto S_state = detectCollision(*sphere1,*sphere2,seconds_type(dt));
            if (S_state.flag == CollisionStateFlag::Collision and (S_state.time>seconds_type(dt_min_Diff) && S_state.time<=dt))
                C.push_back(CollisionObject(sphere1,sphere2,S_state.time));

        }
    }

    for (auto &sphere:_dynamic_spheres){
        for (auto &plane:_static_planes){

            //auto dt_min = sphere->curr_t_in_dt;
            //auto dt_min_Diff = std::max();
            const auto state = detectCollision(*sphere,*plane,dt);

            if (state.flag == CollisionStateFlag::Collision and (state.time>seconds_type(0) && state.time<=dt))
                C.push_back(CollisionObject(sphere,plane,state.time));
        }
    }
        sortAndMakeUnique(C);
        std::reverse(C.begin(),C.end());
        while (!(C.empty())){

            // Take element
            auto c_elem = C.back();
            C.pop_back();

            // Simulate to T
            c_elem.obj1->simulateToTInDt(dt);
            c_elem.obj2->simulateToTInDt(dt);

            // Figure out if objects are static/dynamic
            auto obj1_dsphere = dynamic_cast<DynamicPSphere*>(c_elem.obj1);
//          auto obj1_dcylinder = dynamic_cast<DynamicPSphere*>(c_elem.obj1);
            auto obj2_dsphere = dynamic_cast<DynamicPSphere*>(c_elem.obj2);
            auto obj2_splane = dynamic_cast<const StaticPPlane*>(c_elem.obj2);

            // Compute impact response


             if(obj1_dsphere and obj2_dsphere)
                computeImpactResponse(*obj1_dsphere,*obj2_dsphere,c_elem.t_in_dt);

            else //if(obj1_dsphere and obj2_splane)
                computeImpactResponse(*obj1_dsphere,*obj2_splane,c_elem.t_in_dt);

             // Detect additiona collisions
              if(obj1_dsphere and obj2_dsphere){
                 //std::cout << "Foo"<< std::endl;
                 DynamicColl(obj1_dsphere,obj2_dsphere,c_elem.t_in_dt);
             }
              else if(obj1_dsphere,obj2_splane){
                    std::cout << "Bar"<< std::endl;
                  StaticColl(obj1_dsphere,obj2_splane,c_elem.t_in_dt);
              }
    }
            sortAndMakeUnique(C);
            std::reverse(C.begin(),C.end());

}


void MyController::DynamicColl(DynamicPSphere* d_sphere1,DynamicPSphere* d_sphere2,seconds_type dt)
{
    for(auto &sphere:_dynamic_spheres){

    if(sphere != d_sphere1 and sphere != d_sphere2){
        auto dt_min_Diff = std::max(sphere->curr_t_in_dt,d_sphere1->curr_t_in_dt);

        auto S_state1 = detectCollision(*sphere,*d_sphere1,dt);
        if (S_state1.flag == CollisionStateFlag::Collision and (S_state1.time>seconds_type(dt_min_Diff) && S_state1.time<=dt))
        {
            _collisions.push_back(CollisionObject(sphere,d_sphere1,S_state1.time));
            }

        auto dt_min_Diff2 = std::max(sphere->curr_t_in_dt,d_sphere2->curr_t_in_dt);

        auto S_state2 = detectCollision(*sphere,*d_sphere2,dt);
        if (S_state2.flag == CollisionStateFlag::Collision and (S_state2.time>seconds_type(dt_min_Diff2) && S_state2.time<=dt))
        {
            _collisions.push_back(CollisionObject(sphere,d_sphere2,S_state2.time));
            }
        }
    }
}

void MyController::StaticColl(DynamicPSphere* d_sphere,const StaticPPlane* s_plane,seconds_type dt){

    for(auto &plane:_static_planes){

       if(plane != s_plane){
       auto dt_min = std::max(d_sphere->curr_t_in_dt,seconds_type(0));

        auto state = detectCollision(*d_sphere,*plane,dt);
        if (state.flag == CollisionStateFlag::Collision and (state.time>seconds_type(dt_min) && state.time<=dt))
        {
            _collisions.push_back(CollisionObject(d_sphere,plane,state.time));
            }

      }
    }
}



std::vector<StaticPPlane *> &MyController::getAttachedPlanes(DynamicPSphere *s)
{
    return this->bookingMap[s];
}

void MyController::attachPlaneToSphere(DynamicPSphere* sphere, StaticPPlane* plane)
{
  auto& vector = this->bookingMap[sphere];

  // check if plane is already attached
  auto elem = vector.begin();
  for ( ; *elem != plane; ++elem) {

  }

  // if not, attach it
  if (elem == vector.end())
      vector.push_back(plane);
}



void MyController::localSimulate(double dt)
{
    for(auto &sphere:_dynamic_spheres){
        sphere->curr_t_in_dt = seconds_type(0);
    }
    //add collision algorithm here
     collisionAlgorithm(seconds_type(dt));
    //movement for dynamic spheres

     for(auto &sphere:_dynamic_spheres){

         sphere->simulateToTInDt(seconds_type(dt));
     }
}

std::unique_ptr<Controller> unittestCollisionControllerFactory() {
    return std::make_unique<MyController>();
}

void collision::DynamicPhysObject<GMlib::
PSphere<float> >::setUV(StaticPPlane *_surface, DynamicPSphere *_sphere, float u, float v)
{
    _surface->estimateClpPar(_sphere->getPos(),u,v);
}

GMlib::Vector<float, 3> collision::DynamicPhysObject<GMlib::
PSphere<float> >::getSurfNormal(StaticPPlane *_surface, DynamicPSphere *_sphere, float u, float v)
{
    _surface->getClosestPoint(_sphere->getPos(),u,v);
    GMlib::DMatrix<GMlib::Vector<float,3>> sMatrix = _surface->evaluate(u,v, 1,1);
    GMlib::UnitVector<float,3> norm = sMatrix[0][1] ^ sMatrix[1][0];
    return norm;
}
**/



CollisionState
   detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S0,
                    const DynamicPhysObject<GMlib::PSphere<float>>& S1,
                    seconds_type                                    dt)
   {

       auto max_dt = dt;
       auto min_dt = max(S0.curr_t_in_dt, S1.curr_t_in_dt );
       auto new_dt = max_dt -min_dt;
       auto p0 = S0.getMatrixToScene() *S0.getPos().toType<double>();
       auto r0 = S0.getRadius();
       auto p1 = S1.getMatrixToScene() *S1.getPos().toType<double>();
       auto r1 = S1.getRadius();
       auto r = r0 + r1;
       auto Q = p1 - p0;
       auto R = S1.computeTrajectory(new_dt) - S0.computeTrajectory(new_dt);
       auto a = R * R;
       auto b = Q * R;
       auto c = (Q * Q) - r*r;
       auto epsilon = 1e-6;

       if ((std::abs(c))< epsilon)
       {
           return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityParallelAndTouching));
       }
       else if ((std::abs(a))< epsilon)
       {
           return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityParallel));
       }
       else if ((b*b - a*c)< 0)
       {
           return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityNoCollision));
       }
       else
       {
           auto x = (-b - sqrt(b*b - a*c))/a;
           return (CollisionState (x*new_dt +min_dt,CollisionStateFlag::Collision));
       }

   }

   CollisionState
   detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S,
                    const StaticPhysObject<GMlib::PPlane<float>>&   P,
                    seconds_type                                    dt)
   {

       auto max_dt = dt;
       auto min_dt = S.curr_t_in_dt;
       auto new_dt = max_dt -min_dt;
       auto p = S.getMatrixToScene() * S.getPos();
       auto r = S.getRadius();

       auto unconst_P = const_cast <StaticPhysObject<GMlib::PPlane<float>>&>(P);
       const auto M = unconst_P.evaluateParent(0.5f,0.5f,1,1);
       auto q = M(0)(0);
       auto u = M(1)(0);
       auto v = M(0)(1);
       auto n = GMlib::Vector<float,3>(u ^ v).getNormalized();
       auto d = (q + r * n) - p;

       auto ds = S.computeTrajectory(new_dt);
       auto epsilon = 1e-6;

       if ((std::abs(d * n))< epsilon)
               {
                   return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityParallelAndTouching));
               }

       else if ((std::abs(ds * n))< epsilon)
       {
           return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityParallel));
       }

       else
       {
           auto x = (d * n) / (ds * n);
           return (CollisionState (x*new_dt +min_dt,CollisionStateFlag::Collision));
       }

   }


   CollisionState
   detectCollision (const DynamicPhysObject<GMlib::PSphere<float>>& S,
                    const StaticPhysObject<GMlib::PBezierSurf<float>>& B,
                    seconds_type                                    dt)
   {

       auto max_dt = dt;
       auto min_dt = S.curr_t_in_dt;
       auto new_dt = max_dt -min_dt;
       float u,v,t,delta_u,delta_v,delta_t;
       u = 0.5;
       v = 0.5;
       t = 0.0;
       delta_u = 0.5;
       delta_v = 0.5;
       delta_t = 0.0;
       auto p0 = S.getMatrixToScene() *S.getPos();
       auto r = S.getRadius();
       auto unconst_B = const_cast <StaticPhysObject<GMlib::PBezierSurf<float>>&>(B);
       GMlib::Vector<double, 3> ds = S.computeTrajectory(new_dt);
       GMlib::SqMatrix<double,3> A;

       //iteration

       for (int i=0; i<=5;i++){
       auto p = p0 + ds*t;
       const auto M = unconst_B.evaluateParent(u,v,1,1);
       auto q = M(0)(0);
       const auto Su = M(1)(0);
       const  auto Sv = M(0)(1);
       auto Sn = GMlib::Vector<float,3> (Su ^ Sv).getNormalized();
       A.setCol(Su,0);
       A.setCol(Sv,1);
       A.setCol(-ds,2);
       GMlib::SqMatrix<float,3> A_inv = A;
       A_inv.invert();
       GMlib::APoint<float, 3> b = GMlib::Vector<float, 3> {p-q-Sn*r};
       GMlib::APoint<float, 3> X = A_inv*b;
       delta_u = X(0);
       delta_v = X(1);
       delta_t = X(2);
       u += delta_u;
       v += delta_v;
       t += delta_t;

       auto epsilon = 1e-5;

       if ( (delta_t < epsilon) && (std::abs(delta_u) < epsilon) && (std::abs(delta_v)< epsilon))
       {
           return (CollisionState (seconds_type(delta_t),CollisionStateFlag::Collision));
       }

       }
       return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityNoCollision));

   }


   CollisionState
   detectCollision (const DynamicPSphere&  S,
                    const StaticPCylinder& C,
                    seconds_type dt)
{

   auto max_dt = dt;
   auto min_dt = S.curr_t_in_dt;
   auto new_dt = max_dt -min_dt;
   auto p = S.getMatrixToScene() * S.getPos();
   auto r_s = S.getRadius();
   auto r_c = S.getRadius();
   auto unconst_C = const_cast <StaticPhysObject<GMlib::PCylinder<float>>&>(C);
   const auto M = unconst_C.evaluateParent(0.5f,0.5f,1,1);
   auto q = M(0)(0);
   auto u = M(1)(0);
   auto v = M(0)(1);
   auto n = GMlib::Vector<float,3>(u ^ v).getNormalized();
   auto d = (q + (r_s + r_c) * n) - p;

   auto ds = S.computeTrajectory(new_dt);
   auto epsilon = 1e-6;

   if ((std::abs(d * n))< epsilon)
   {
      return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityParallelAndTouching));
   }

   else if ((std::abs(ds * n))< epsilon)
   {
       return(CollisionState(seconds_type (0.0),CollisionStateFlag::SingularityParallel));
   }

   else
   {
       auto x = (d * n) / (ds * n);
       return (CollisionState (x*new_dt +min_dt,CollisionStateFlag::Collision));
   }

}



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
       auto S1_new_vel = v1_n * n + new_v1_d * normal_d;
       S0.velocity = S0_new_vel;
       S1.velocity = S1_new_vel;

   }

   void
   computeImpactResponse (DynamicPhysObject<GMlib::PSphere<float>>& S,
                          const StaticPhysObject<GMlib::PPlane<float>>&   P,
                          seconds_type                              dt)
   {

//       auto unconst_P = const_cast <StaticPhysObject<GMlib::PPlane<float>>&>(P);
//       const auto M = unconst_P.evaluateParent(0.5f,0.5f,1,1);
//       auto u = M(1)(0);
//       auto v = M(0)(1);
//       auto n = GMlib::Vector<float,3>(u ^ v).getNormalized();
//       auto vel = S.velocity * n;

//       auto new_velocity = S.velocity - ((2* vel)*n)*0.90;
//       S.velocity = new_velocity ;

       auto &unconst_P = const_cast<StaticPhysObject<GMlib::PPlane<float>>&>(P);
       auto s_pos =S.getPos();
       auto _Radius = S.getRadius();
       auto plane_pos = unconst_P.evaluateParent(0.5f, 0.5f,1, 1);

       auto v = plane_pos(0)(1);
       auto u = plane_pos(1)(0);
       auto n = u^v;
       auto nNormal = GMlib::Vector<float,3>(n).getNormalized();

       auto newVel = (S.velocity - 2*(S.velocity*nNormal)*nNormal)*0.90;
       S.velocity = newVel;

   }


   std::unique_ptr<Controller> unittestCollisionControllerFactory(){ return std::make_unique<MyController> (); }

   GMlib::Vector<double,3> DynamicPhysObject<GMlib::PSphere<float> >::externalForces() const {
       assert (environment != nullptr);
       return this->environment->externalForces().toType<double>();

   }


   void MyController::localSimulate(double dt) {

       //we need to reset currentTInDt of all dynamic objects

       for(auto& sphere : _dynamic_spheres){
           sphere->curr_t_in_dt =seconds_type(0.0);
       }
       //state changes
       //detectStateChanges(dt);

/**
       for(auto& sphere : _dynamic_spheres){
          for(auto& plane: _static_planes){

            //detect singularities (Collisions, Parallel and Touching, Parallel, NoCollision)
            auto col = detectCollision(*sphere,*plane, seconds_type(dt));

            if (col.flag == CollisionStateFlag::SingularityParallelAndTouching &&
                    !( std::find(attachedPlanes[sphere].begin(), attachedPlanes[sphere].end(), plane) != attachedPlanes[sphere].end() )){
                std::cout << "Attached" << std::endl;
                setAttachedObjects(sphere,plane);
                sphere->addPlanes(plane);}

//            if (col.flag != CollisionStateFlag::SingularityParallelAndTouching &&
//            ( std::find(attachedPlanes[sphere].begin(), attachedPlanes[sphere].end(), plane) != attachedPlanes[sphere].end() ))
//                    // !(attachedPlanes[sphere].empty())){
//                //attachedPlanes[sphere].erase(plane);
//            }
            for (auto sphere = _dynamic_spheres.begin() ; sphere != _dynamic_spheres.end() ; ++sphere){
                    stateChangeObject stateObject = detectStateChanges(*sphere,dt);
                    (*sphere)->state= stateObject.stateChanges;
                }

            //if no collisions
            if (col.flag == CollisionStateFlag::SingularityNoCollision){
                     attachedPlanes[sphere].erase(plane);}

            //if parallel
 //           if (col.flag == CollisionStateFlag::SingularityParallel){
              //if distance from the sphere to plane > radius of sphere,
                //sphere is free and flying
//            }


          }
       }
**/

//singularities
/**
//       detectSingularities(dt);
//       std::reverse(_singularities.begin(), _singularities.end());

//       for( auto& singularity : _singularities) {

//           detectStateChange(singularity);
//       }
**/

    //detect statechanges old
/**
       for (auto iter = _dynamic_spheres.begin(); iter != _dynamic_spheres.end(); iter++){

            auto stateObject = detectStateChanges(*iter,dt);

            (*iter)->state= stateObject.stateChanges;

        }
**/
       //Detect collisions
       detectCollisions(dt);
       sortAndMakeUnique(_collisions);

//&& !_stateChanges.empty())
       while (!_collisions.empty()){ //&& !_singularities.empty()){

           auto col = _collisions.begin();
           auto col_time =  col->first;

//singularities
//           auto singularity = _singularities.begin();
 //          auto sing_time =  singularity->first;

//collision happens first
   //if (col_time<sing_time){

           auto col_obj1 = col->second.obj1;
           auto col_obj2 = col->second.obj2;

           _collisions.erase(col); //col is now invalid but it won't be used until next iteration when we redefine it
           if (auto sphere2 = dynamic_cast<DynamicPSphere*> (col_obj2)){
               auto sphere1 = dynamic_cast<DynamicPSphere*> (col_obj1);
               sphere1->simulateToTInDt(col_time);
               sphere2->simulateToTInDt(col_time);
               computeImpactResponse(*sphere1,*sphere2,col_time);
           }

           else {
               auto sphere = dynamic_cast<DynamicPSphere*> (col_obj1);
               auto plane = dynamic_cast<StaticPPlane*> (col_obj2);
               sphere->simulateToTInDt(col_time);
               computeImpactResponse(*sphere,*plane,col_time);
           }

           //Detect more collisions
           detectCollisions(dt);
           sortAndMakeUnique(_collisions);

           //detectStateChanges(dt);//detect state changes
            //}//if ends



   //state changes happens first
//            else{
//                singularity->second.sphere_obj->simulateToTInDt(sing_time);
//                singularity->second.sphere_obj->state = singularity->second.stateChange;
//                _singularities.erase(singularity);
//            }

   //Detect more collisions
//            detectCollisions(dt);
//            sortAndMakeUnique(_collisions);

//              //Detect more state changes
//            detectStateChanges(dt);

          }
       for(auto& sphere : _dynamic_spheres){
           sphere->simulateToTInDt(seconds_type(dt));
       }

   }
   void DynamicPhysObject<GMlib::PSphere<float> >::simulateToTInDt( seconds_type t) {

       auto t0 = seconds_type(t - this->curr_t_in_dt);
       auto Mi = this->getMatrixToSceneInverse();
       //move
       //update a ds
       auto r = this->getRadius();
       auto s = this->getCenterPos();

       auto ds = this->computeTrajectory(t0);
       auto p = s+ds;
       GMlib::Vector <float,3>n {0.0f,0.0f,0.0f};

       for (auto &it :_planes){
            //std::cout << _planes.size();
             auto M = it->evaluateParent(0.5f,0.5f,1,1);
             auto u = M(1)(0);
             auto v = M(0)(1);
             auto n = GMlib::Vector<float,3>(u ^ v);
             auto normal = GMlib::Vector<float,3>(u ^ v);
              n+=normal;
        }
       n= GMlib::Vector <float,3>(n/_planes.size()).getNormalized();
       if (this->state == states::Rolling) {
           // ...
           ds = ds+r+(p*s)*n;

       }
//       else
//         ds = this->computeTrajectory(t0);

       this->translateParent(Mi*ds);
       this->curr_t_in_dt =t;
       //update physics
       auto F = this->externalForces();
       auto c = t0.count();
       auto a = F*c;
       this->velocity += a;

   }

   void MyController::detectCollisions(double dt){

       //loop for collision between dynamic objects (only spheres for now)
       for (auto it1 = _dynamic_spheres.begin() ; it1 != _dynamic_spheres.end() ; ++it1){
           for (auto it2 = it1+1 ; it2 != _dynamic_spheres.end() ; ++it2 ){
               auto col = collision::detectCollision(**it1,**it2,seconds_type(dt));

               const auto &sphere1= *it1;
               const auto &sphere2= *it2;

               auto min_ctidt = std::max(sphere1->curr_t_in_dt, sphere2->curr_t_in_dt);

               if (col.flag == CollisionStateFlag::Collision && col.time < seconds_type(dt) && col.time > min_ctidt ){
                   _collisions.emplace(col.time,CollisionObject(sphere1,sphere2,col.time));
               }
           }
       }

       //loop for collision with static objects (only dynamic spheres with static planes for now)
       //add setAttached for flag Parallel and toching
       for (auto &it1 : _dynamic_spheres){
           for (auto &it2 : _static_planes){

               auto col = collision::detectCollision(*it1,*it2,seconds_type(dt));

               if (col.flag == CollisionStateFlag::Collision && col.time < seconds_type(dt) && col.time > it1->curr_t_in_dt){
                   _collisions.emplace(col.time,CollisionObject(it1,it2,col.time));}
               }
           }
       }

   std::vector<StaticPPlane*> const MyController::getAttachedPlanes(DynamicPSphere* sphere) {

       return (attachedPlanes[sphere]);

   }


   void MyController::detectStateChanges(double dt)
   {
       for (auto iter = _dynamic_spheres.begin() ; iter != _dynamic_spheres.end() ; ++iter){
           auto singularity = detectStateChange(*iter,dt);
           const auto &sphere= *iter;
           if (singularity.stateChange != states::NoChange){
               _singularities.emplace(sphere->curr_t_in_dt,stateChangeObj(sphere,_static_planes,sphere->curr_t_in_dt, singularity.stateChange));
           }
       }

   }

   void MyController::detectSingularities(double dt)
   {
/**
//    for( auto& sphere : _dynamic_spheres) {
//         for( auto& plane : _static_planes) {

//         auto r = sphere->getRadius();
//         auto p = sphere->getMatrixToScene() * sphere->getPos();

//         auto epsilon = 1e-5;

//         auto dts = seconds_type(dt);
//         auto max_dt = dts;
//         auto min_dt = sphere->curr_t_in_dt;
//         auto new_dt = max_dt -min_dt;
//         auto ds = sphere->computeTrajectory(new_dt);

//         auto M = plane->evaluateParent(0.5f,0.5f,1,1);
//         auto q = M(0)(0);
//         auto u = M(1)(0);
//         auto v = M(0)(1);
//         auto n = GMlib::Vector<float,3>(u ^ v);

//                      // Check for formulas
//         auto d = (q + r * n) - p;
//         auto still = std::abs(((-n*r) * ds) -(ds*ds));
//         auto dsn = ds * n;
//         auto dn = d*n;

//            if ( std::abs( ((-n*r) * ds) -(ds*ds) ) < epsilon ) {

//                auto singularity = stateChangeObj(sphere, plane, seconds_type(dt), states::Still);
//                _singularities.push_back(singularity);
//                 }
//             else if(std::abs(ds * n) > 0 && std::abs(((-n*r) * ds) -(ds*ds)) > epsilon) { //&& !(attachedPlanes[sphere].empty())

//                 auto singularity = stateChangeObj(sphere, plane, seconds_type(dt), states::Rolling);
//                 _singularities.push_back(singularity);
//                  }
//              else if( std::abs(ds * n) <= 0) {

//                 auto singularity = stateChangeObj(sphere, plane, seconds_type(dt), states::Free);
//                 _singularities.push_back(singularity);
//                      }
//                  }
//              }
**/
   }

   stateChangeObj MyController::detectStateChange(DynamicPSphere *sphere, double dt)
   {
       std::vector<StaticPPlane*> p = _static_planes;
       auto r = sphere->getRadius();
       auto pos = sphere->getMatrixToScene() * sphere->getPos();

       auto epsilon = 1e-5;
       auto max_dt = seconds_type(dt);
       auto min_dt = sphere->curr_t_in_dt;

       auto new_dt = max_dt -min_dt;
       auto ds = sphere->computeTrajectory(new_dt);
       auto attached_planes = attachedPlanes[sphere] ;

       GMlib::APoint<float,3> q;
       GMlib::Vector <float,3> n {0.0f,0.0f,0.0f};

       if (attached_planes.empty()){// sphere is free
           for (auto &plane : _static_planes){
               auto M = plane->evaluateParent(0.5f,0.5f,1,1);
               auto q= M(0)(0);
               auto u = M(1)(0);
               auto v = M(0)(1);
               auto n = GMlib::Vector<float,3>(u ^ v).getNormalized();
               auto d = (q + r * n) - pos;

               if (std::abs(((-n*r) * ds) -(ds*ds)) < epsilon){
                    p.push_back(plane);
                               // return stateChangeObject(sphere,states::Still);
                    return stateChangeObj(sphere, p,min_dt,states::Still);
                    }
               else if ( std::abs(d*n)< epsilon && ds * n <= 0 ){
                    p.push_back(plane);
                               // return stateChangeObject(sphere,  states::Rolling);
                    return stateChangeObj(sphere, p,min_dt,states::Rolling);
                   }
               else stateChangeObj(sphere, p,min_dt,states::NoChange);
          }
       }

       else{ //attached
              for (auto &plane :attached_planes){
                   auto M = plane->evaluateParent(0.5f,0.5f,1,1);
                   auto pos= M(0)(0);
                   auto u = M(1)(0);
                   auto v = M(0)(1);
                   auto normal = GMlib::Vector<float,3>(u ^ v);
                   n+=normal;
                   q=pos;
               }
               n= GMlib::Vector <float,3>(n/attached_planes.size()).getNormalized();

               auto d = (q + r * n) - pos;
               auto bla=std::abs(((-n*r) * ds) -(ds*ds));
               auto dsn= ds * n;
               auto dn= d*n;

               if (sphere->state == states::Rolling){
                    if (ds * n > 0){
                        std::cout<<"state changes from Rolling to Free";
                           //detachement
                        attachedPlanes[sphere].clear(); //empty vector
                           //return (stateChangeObject(sphere, states::Free));
                        return stateChangeObj(sphere, attachedPlanes[sphere], min_dt, states::Free);
                    }
                    else if (std::abs(((-n*r) * ds) -(ds*ds)) < epsilon){
                        std::cout<<"state changes from Rolling to Still";
                           //return (stateChangeObject(sphere,states::Still)) ;
                        return stateChangeObj(sphere, p,min_dt,states::Still);
                       }
                    else stateChangeObj(sphere, p,min_dt,states::NoChange);
                   }

               else if (sphere->state == states::Still){
                    if (std::abs(((-n*r) * ds) -(ds*ds)) > epsilon){
                           //Correct trajectory
                       std::cout<<"state changes from still to Rolling";
                           //return (stateChangeObject(sphere,  states::Rolling));
                       return stateChangeObj(sphere, p,min_dt,states::Rolling);
                       }
                    else if (ds * n > 0  ){
                       std::cout<<"state changes from still to Free";
                           //detachement
                           attachedPlanes[sphere].clear(); //empty vector
                          // return (stateChangeObject(sphere,  states::Free));
                       return stateChangeObj(sphere, attachedPlanes[sphere],min_dt, states::Free);
                       }
                    else stateChangeObj(sphere, p,min_dt,states::NoChange);
                   }
                }


   }

   void MyController::setAttachedObjects(DynamicPSphere *sphere, StaticPPlane *plane){

//       for (auto iter = attachedPlanes[sphere].begin() ; iter != attachedPlanes[sphere].end() ; ++iter)
//       //const bool is_in = attachedPlanes[sphere].find(plane) != attachedPlanes[sphere].end();
//      if (!(is_in))
       attachedPlanes[sphere].push_back(plane);
       std::sort(attachedPlanes[sphere].begin(),attachedPlanes[sphere].end());
       std::unique(attachedPlanes[sphere].begin(),attachedPlanes[sphere].end());


   }

//make deattachment
   void MyController::setDeattachedObjects(DynamicPSphere *sphere, StaticPPlane *plane)
   {
       std::swap(plane,attachedPlanes[sphere].back());
       attachedPlanes[sphere].pop_back();

   }

   void setState(DynamicPSphere* sphere, states state){

       sphere->state =state;

   }

   states getState(DynamicPSphere* sphere){

       return sphere->state;

   }


   void MyController::add (DynamicPSphere* const sphere) {
       _dynamic_spheres.push_back(sphere);
       auto m=0;
       StaticPPlane* p;
       //setAttachedObjects(sphere,_static_planes.back());
       sphere->environment = &_environment;
   }

   void collision::DynamicPhysObject<GMlib::PSphere<float> >::addPlanes(StaticPPlane *plane)
   {
       this->_planes.insert(plane);
   }

   GMlib::Vector<double, 3> collision::DynamicPhysObject<GMlib::PSphere<float> >::adjustedTrajectory(seconds_type dt)
   {
       auto r = this->getRadius();
       auto s = this->getMatrixToScene() * this->getPos();
       auto ds = this->computeTrajectory(dt);
       auto p = s+ds;
       //dont forget to add a planes
       auto planes = _planes;
       GMlib::Vector <float,3>n {0.0f,0.0f,0.0f};

       for (auto &it :planes){
           auto M = it->evaluateParent(0.5f,0.5f,1,1);
           auto q = M(0)(0);
           auto u = M(1)(0);
           auto v = M(0)(1);
           auto normal = GMlib::Vector<float,3>(u ^ v);
           n+=normal;
      }
       n= GMlib::Vector <float,3>(n/planes.size()).getNormalized();

      auto dsAdjusted = ds - ds*n*n;

      return dsAdjusted;
   }








} // END namespace collision

