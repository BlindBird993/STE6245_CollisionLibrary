#ifndef COLLISION_LIBRARY_H
#define COLLISION_LIBRARY_H

// collision library interface
#include <collision_interface.h>
#include <cmath>
#include <unordered_map>


namespace collision
{
/**
//states
enum class States{
    NoChange,
    Free,
    Rolling,
    Still
};

struct StateChangeObject {
    DynamicPSphere* obj1;
    States stateType;

    StateChangeObject(DynamicPSphere* obj1,States stateType): obj1{obj1},stateType{States::Free} {}

};


    template <>
    class DynamicPhysObject<GMlib::PSphere<float>> : public DynamicPhysObject_Base<GMlib::PSphere<float>> {
    public:
        using DynamicPhysObject_Base<GMlib::PSphere<float>>::DynamicPhysObject_Base;

        States                                                                     state= States::Free;

        void    simulateToTInDt( seconds_type t );
        void    NewSimulateToTInDt( seconds_type t,StaticPPlane* _surface);
        GMlib::Vector<double, 3> externalForces () const;
        void computeStep(seconds_type dt);

        // make attached
        void setUV(StaticPPlane* _surface,DynamicPSphere* _sphere,float u, float v);
        GMlib::Vector<float,3> getSurfNormal(StaticPPlane* _surface,DynamicPSphere* _sphere,float u, float v);

        //variable for storing ds


        void setAttached(StaticPPlane* _p);
        StaticPPlane* getAttachedTo();
        StaticPPlane* getPlane();
        void          setPlane(StaticPPlane* _surf);

        const Controller* _controller;

        void correctTrajectory(GMlib::Vector<double,3> updated_ds);

        struct CollisionState {
            seconds_type       time;
            CollisionStateFlag flag;
            CollisionState( seconds_type t, CollisionStateFlag f = CollisionStateFlag::Collision )
                : time{std::move(t)}, flag{std::move(f)} {}
        };

    public:
        //States          state{States::Free};
        StaticPPlane*  _plane;

        //new ds for attach surfaces
        mutable GMlib::Vector<double,3>  new_ds;

    };

    class MyController : public Controller {
        GM_SCENEOBJECT(MyController)
        public:
//new


        explicit MyController () = default;

        void detectCollisions(double dt);


//            enum class States{
//                Free,
//                Sliding,
//                AtRest,
//                NoChange
//            };


//            struct StateChangeObject {
//                DynamicPSphere* obj1;
//                States stateType;

//                StateChangeObject(DynamicPSphere* obj1,States stateType): obj1{obj1},stateType{States::Free} {}

//            };

//MyController add sphere
           void add (DynamicPSphere* const sphere); //{
//                    sphere->environment = &_environment;
//                    _dynamic_spheres.push_back(sphere);
//           }
            void add (StaticPSphere* const sphere) { _static_spheres.push_back(sphere); }
            void add (StaticPPlane* const plane) { _static_planes.push_back(plane); }
            void add (StaticPCylinder* const cylinder) { _static_cylinders.push_back(cylinder); }
            void add (StaticPBezierSurf* const surf) { _static_bezier_surf.push_back(surf); }


            void collisionAlgorithm (seconds_type dt);
            void DynamicColl (DynamicPSphere* d_sphere1,DynamicPSphere* d_sphere2,seconds_type dt);
            void StaticColl (DynamicPSphere* d_sphere,const StaticPPlane* plane, seconds_type dt);
            void detectCol(DynamicPSphere* sphere, seconds_type dt);

            std::vector<StaticPPlane*> &getAttachedPlanes(DynamicPSphere* s);

            std::unordered_map<DynamicPSphere*,std::vector<StaticPPlane*>> bookingMap;

            void attachPlaneToSphere(DynamicPSphere* sphere, StaticPPlane* plane);

            StateChangeObject detectStateChanges( DynamicPSphere* sphere,double dt);

            States getState(DynamicPSphere* sphere);

            void setState(DynamicPSphere* sphere, States state);

            void setAttachedObjects(DynamicPSphere*  sphere  , StaticPPlane* plane);



        protected:
            std::vector<DynamicPSphere*>    _dynamic_spheres;

            std::vector<StaticPSphere*>     _static_spheres;

            std::vector<StaticPPlane*>      _static_planes;

            std::vector<StaticPCylinder*>   _static_cylinders;
            std::vector<StaticPBezierSurf*> _static_bezier_surf;
            collision::DefaultEnvironment _environment;

            std::multimap<seconds_type,collision::CollisionObject> _collisions;

            //std::vector<collision::CollisionObject> _collisions;

            void localSimulate(double dt) override final;
        };

    template <class PSurf_T, typename... Arguments>
    std::unique_ptr<DynamicPhysObject<PSurf_T>> unittestDynamicPhysObjectFactory(Arguments... parameters) {

        return std::make_unique<DynamicPhysObject<PSurf_T>>(parameters...);
    }

    template <class PSurf_T, typename... Arguments>
    std::unique_ptr<StaticPhysObject<PSurf_T>> unittestStaticPhysObjectFactory(Arguments... parameters) {

     return std::make_unique<StaticPhysObject<PSurf_T>>(parameters...);
    }


    template <class Container_T >
    void sortAndMakeUnique(Container_T &container)
    {
        //sort
        std::sort(container.begin(),container.end(),[](const auto &a, const auto &b){
            return a.t_in_dt < b.t_in_dt;
        });

        //make uique
        auto pred = [](const auto& a, const auto& b){

            auto is_dynamic = [](const auto *obj){
                if (dynamic_cast<const DynamicPSphere*>(obj)) return true;

                return false;
          };
          if(a.obj1 == b.obj2) return true;
          if (a.obj2 == b.obj1) return true;
          if (a.obj1 == b.obj1) return true;
          if ( ( is_dynamic(a.obj2) or is_dynamic(b.obj2) )
               and a.obj2 == b.obj2 ) return true;

       return false;
   };
        typename Container_T::iterator NE=std::end(container);
        for (auto itr = std::begin(container);itr!=NE;++itr){

            for (auto ritr = NE-1;ritr != itr;--ritr){

                if ((pred(*itr,*ritr))){

                    std::swap(*ritr,*(NE-1));
                    NE--;
                }
            }

    }
     container.erase(NE,std::end(container));
 }



    // Sort
    {
      std::sort( begin(container), end(container), [](const auto& a, const auto& b) {
          return a.t_in_dt < b.t_in_dt;
      });

      // Make unique

      auto pred =  [](const auto &a, const auto &b) {

          auto is_d_pred = []( const auto* obj ) {
            if(dynamic_cast<const DynamicPSphere*>(obj)) return true;



            return false;
          };

          if( a.obj1 == b.obj1 ) return true;
          if( a.obj1 == b.obj2 ) return true;
          if( a.obj2 == b.obj1 ) return true;
          if( ( is_d_pred(a.obj2) or is_d_pred(b.obj2) )
                  and a.obj2 == b.obj2 ) return true;

          return false;
      };

      typename Container_T::iterator NE = std::end(container);
      for( auto first_iter = std::begin(container); first_iter != NE; ++first_iter) {

          for( auto r_iterator = NE - 1; r_iterator != first_iter; --r_iterator) {

                  if( (pred(*first_iter, *r_iterator))) {
                          std::swap( *r_iterator, *(NE-1) );
                          NE--;

                      }
              }
          }

          container.erase( NE, std::end( container ) );

//take one, compares to all the others, mark for removal and again...
//    template <class PSurf_T, typename... Arguments>
//    std::unique_ptr<StaticPhysObject<PSurf_T>> unittestStaticPhysObjectFactory(Arguments... parameters) {

//        return std::make_unique<StaticPhysObject<PSurf_T>> (parameters...);
//    }


**/




enum class states {
      NoChange,
      Free,
      Rolling,
      Still
  };

//  struct stateChangeObject{
//       DynamicPSphere* obj;
//       states stateChanges;


//      stateChangeObject( DynamicPSphere* object, states s  )
//          : obj{object}, stateChanges{s} {}
//  };

  struct stateChangeObj {
      DynamicPSphere*                 sphere_obj;   // Object whos state will change
      std::vector<StaticPPlane*>      planes;   // Object that obj1 will change state ACCORDING to
      seconds_type                    time;   // Time of singularity
      states                          stateChange;  // State that obj1 will change to

      stateChangeObj
      (DynamicPSphere* obj1, std::vector<StaticPPlane*> obj2, seconds_type t,states s) :
          sphere_obj{obj1}, planes{obj2}, time{t} ,stateChange{s} {}
  };


  class MyController : public Controller {
      GM_SCENEOBJECT (MyController)

  public:

      //explicit MyController () = default;

      void add (DynamicPSphere* const sphere);
      void add (StaticPSphere* const sphere) { _static_spheres.push_back(sphere);  }
      void add (StaticPPlane* const plane) { _static_planes.push_back(plane); }
      void add (StaticPCylinder* const cylinder) { _static_cylinders.push_back(cylinder); }
      void add (StaticPBezierSurf* const surf) { _static_bezier_surf.push_back(surf); }

  protected:
      void localSimulate (double dt) override;
      void detectCollisions(double dt);

      void setAttachedObjects(DynamicPSphere*  sphere  , StaticPPlane* plane);
      void setDeattachedObjects(DynamicPSphere*  sphere  , StaticPPlane* plane);



      std::vector<StaticPPlane*>  const getAttachedPlanes(DynamicPSphere* sphere) ;

      void setState(DynamicPSphere* sphere, states state);

      states getState(DynamicPSphere* sphere);

      void detectStateChanges(double dt);
      stateChangeObj detectStateChange( DynamicPSphere* sphere,double dt);

      void attachPlane(DynamicPSphere*  sphere  , StaticPPlane* plane);
      void detachPlane(DynamicPSphere*  sphere  , StaticPPlane* plane);

      void detectSingularities(double dt);



      std::vector<DynamicPSphere*>                                                                          _dynamic_spheres;
      std::vector<StaticPSphere*>                                                                           _static_spheres;
      std::vector<StaticPPlane*>                                                                            _static_planes;
      std::vector<StaticPCylinder*>                                                                         _static_cylinders;
      std::vector<StaticPBezierSurf*>                                                                       _static_bezier_surf;

      std::multimap<seconds_type,collision::CollisionObject>                                                _collisions;
      std::multimap<seconds_type,stateChangeObj>                                                            _singularities;

      std::unordered_map<DynamicPSphere* , std::vector<StaticPPlane*>>                                      attachedPlanes;
      DefaultEnvironment                                                                                    _environment;

  };

  template <>
  class DynamicPhysObject<GMlib::PSphere<float>> : public DynamicPhysObject_Base<GMlib::PSphere<float>> {
  public:
      using DynamicPhysObject_Base<GMlib::PSphere<float>>::DynamicPhysObject_Base;

      states                                                                     state= states::Free;
      GMlib::Vector<float,3>                                                     trajectory;
      std::set<StaticPPlane*>                                                    _planes;
      void addPlanes(StaticPPlane* plane);
      GMlib::Vector<double,3> adjustedTrajectory(seconds_type dt);


      void    simulateToTInDt( seconds_type t) override;

      GMlib::Vector<double,3> computeTrajectory( seconds_type dt) const override {

              auto t=dt.count();
              auto F = this->externalForces();
              const auto a = 0.5*F*std::pow(t,2)*this->mass;
              return this->velocity*t + a;
      }


      GMlib::Vector<double, 3> externalForces () const override ;

  };

  template <class PSurf_T, typename... Arguments>
  std::unique_ptr<DynamicPhysObject<PSurf_T>> unittestDynamicPhysObjectFactory(Arguments... parameters) {

      return std::make_unique<DynamicPhysObject<PSurf_T>> (parameters...);

  }

  template <class PSurf_T, typename... Arguments>
  std::unique_ptr<StaticPhysObject<PSurf_T>> unittestStaticPhysObjectFactory(Arguments... parameters) {

      return std::make_unique<StaticPhysObject<PSurf_T>> (parameters...);
  }

  template <class Container_T >
  void sortAndMakeUnique( Container_T& container) {

  //make unique;
      for (auto it1 = container.begin() ; it1 != container.end() ; ++it1){

          for (auto it2 = std::next(it1,1) ; it2 != container.end() ;){

              if (it2->second.obj1 == it1->second.obj1 || it2->second.obj2 == it1->second.obj1 ||
                      it2->second.obj1 == it1->second.obj2 ||
                      (it2->second.obj2 == it1->second.obj2 && dynamic_cast<DynamicPSphere*> (it1->second.obj2)))
                    {
                      container.erase(it2++);
                    }
              else
                    {
                      ++it2;
                    }
          }
      }



  }


  template <class Container_T >
  void statesSortAndMakeUnique(Container_T &container)
  {
      //sort
      std::sort(container.begin(),container.end(),[](const auto &a, const auto &b){
          return a.t_in_dt < b.t_in_dt;
      });

      //make uique
      auto pred = [](const auto& a, const auto& b){

          auto is_dynamic = [](const auto *obj){
              if (dynamic_cast<const DynamicPSphere*>(obj)) return true;

              return false;
        };
        if(a.obj1 == b.obj2) return true;
        if (a.obj2 == b.obj1) return true;
        if (a.obj1 == b.obj1) return true;
        if ( ( is_dynamic(a.obj2) or is_dynamic(b.obj2) )
             and a.obj2 == b.obj2 ) return true;

     return false;
 };
      typename Container_T::iterator NE=std::end(container);
      for (auto itr = std::begin(container);itr!=NE;++itr){

          for (auto ritr = NE-1;ritr != itr;--ritr){

              if ((pred(*itr,*ritr))){

                  std::swap(*ritr,*(NE-1));
                  NE--;
              }
          }

  }
   container.erase(NE,std::end(container));
}














} // END namespace collision



#endif //COLLISION_LIBRARY_H
