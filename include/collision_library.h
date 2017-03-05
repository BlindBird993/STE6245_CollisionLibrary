#ifndef COLLISION_LIBRARY_H
#define COLLISION_LIBRARY_H

// collision library interface
#include <collision_interface.h>

#include <typeinfo>
#include <unordered_map>


namespace collision
{
#define States_enum_class_and_StateChangeObj {
//States
enum class States {
    Free,
    Rolling,
    Still,
};


//struct StateChangeObj;
struct StateChangeObj {
    DynamicPSphere*                             obj;         // Object whos state will change
    std::unordered_set<StaticPPlane*>           planes;     // Object that obj1 will change state ACCORDING to
    seconds_type                                time;               // Time of singularity
    States                                      state;              // State that obj1 will change to

    StateChangeObj
    (DynamicPSphere* o1, std::unordered_set<StaticPPlane*> planes, seconds_type t, States s) :
        obj{o1}, planes{planes}, time{t} ,state{s} {}
};
#define end1 }

#define MyControllerClass {
class MyController : public Controller {
    GM_SCENEOBJECT (MyController)
    public:

        explicit MyController () = default;

    void add (DynamicPSphere* const sphere);
    void add (StaticPSphere* const sphere);
    void add (StaticPPlane* const plane);
    void add (StaticPCylinder* const cylinder);
    void add (StaticPBezierSurf* const surf);

    // Change "value" and return type if the sphere should be able to be attached to objects other than planes
    std::unordered_set<StaticPPlane *> getAttachedObjects(DynamicPSphere* sphere);


    void setAttachedObjects( std::unordered_set<StaticPPlane*> planes, DynamicPSphere* sphere );
    void detachObjects(DynamicPSphere* sphere);

    void reverseMethod(std::vector<StateChangeObj> singularitiesContainer,std::vector<collision::CollisionObject> collisionsContainer);

    void
    detectStateChanges(double dt);//fill the container of singularities

    StateChangeObj
    detectStateChange(DynamicPSphere* sphere, double dt); //check the state

    void
    detectStates (StateChangeObj& state, double dt);//update the state

    GMlib::Vector<float,3>
    closestPoint(DynamicPSphere* sphere, seconds_type dt);//get closest point

    void detectCollisions ( collision::CollisionObject& col, double dt);//checks what objects are collided

    Environment                                 noGravEnv;//no gravity environment
    DefaultEnvironment                          env;// default environment

    std::unordered_map<DynamicPSphere*, std::unordered_set<StaticPPlane*>>    _attachedPlanes;//container for attached planes



private:

    template <typename T_s>
    void sphereCollision(T_s* sphere,seconds_type dt);//collison detection

    template <class Container>
    void sortAndMakeUniqueStates(Container& c);//sort and make states

    template <typename Container_1, typename Container_2>
    void crossUnique( Container_1 container1, Container_2 container2);// take two sorted and unique containers and crosscheck them,
    //update containers

protected:
    void localSimulate (double dt) override;//simulation loop
//containers
    std::vector<DynamicPSphere*>                _dynamic_spheres;
    std::vector<StaticPSphere*>                 _static_spheres;
    std::vector<StaticPPlane*>                  _static_planes;
    std::vector<StaticPCylinder*>               _static_cylinders;
    std::vector<StaticPBezierSurf*>             _static_bezier_surf;

    std::vector<collision::CollisionObject>     _collisions;
    std::vector<StateChangeObj>                 _singularities;
};

#define end1 }

#define DynamicPhysObjectSphereClass {
template <>
class DynamicPhysObject<GMlib::PSphere<float>> : public DynamicPhysObject_Base<GMlib::PSphere<float>> {
public:
    using DynamicPhysObject_Base<GMlib::PSphere<float>>::DynamicPhysObject_Base;


    MyController*   _sphereController; //point to controller
    States          _state = States::Free;     // default state

    //methods to change position of the sphere
    void moveUp();
    void moveDown();
    void moveLeft();
    void moveRight();


    void translateUp();
    void translateDown();
    void translateLeft();
    void translateRight();

    //get and set velocity
    void setVelocity(const GMlib::Vector<double,3> vel);
    GMlib::Vector<double,3> getVelocity();

    GMlib::Vector<double,3>
    adjustedTrajectory (seconds_type dt);// adjust trajectory


    void            simulateToTInDt( seconds_type t ) override;

    GMlib::Vector<double, 3>
    computeTrajectory (seconds_type dt) const override; //compute trajectory

    GMlib::Vector<double, 3> externalForces () const override;

};
#define end1 }

#define PlaneTemplate {
template <>
class StaticPhysObject<GMlib::PPlane<float>> :  public PhysObject<GMlib::PPlane<float>, PhysObjectType::Static> {
public:
    int id=0;

    using PhysObject<GMlib::PPlane<float>, PhysObjectType::Static>::PhysObject;


    void simulateToTInDt (seconds_type) override {}

    int getId() const;
    void setId(int value);


};
#define end1 }

template <class PSurf_T, typename... Arguments>
std::unique_ptr<DynamicPhysObject<PSurf_T>> unittestDynamicPhysObjectFactory(Arguments... parameters) {

    return std::make_unique<DynamicPhysObject<PSurf_T>> (parameters...);

}

template <class PSurf_T, typename... Arguments>
std::unique_ptr<StaticPhysObject<PSurf_T>> unittestStaticPhysObjectFactory(Arguments... parameters) {

    return std::make_unique<StaticPhysObject<PSurf_T>> (parameters...);
}

#define crossUniqueMethod {
/**
  *
  * Code was written with help from Bjørn-Richard Pedersen
  **/
// crosscheck for states and collisions
template <typename Container_1, typename Container_2>
void MyController::crossUnique(Container_1 container1, Container_2 container2) {

    std::vector<collision::CollisionObject>     realColContainer; //new containers
    std::vector<StateChangeObj>                 realSingContainer;

    auto objPred =  [](const auto &a, const auto &b) {

        if( a.obj1 == b.obj or a.obj2 == b.obj ) return true; //check if equal

        return false;
    };

    auto timePred =  [](const auto &a, const auto &b) {

        if( a.t_in_dt < b.time ) return true; //compare objects time, if collision is before the state or not

        return false;//if state is first
    };

    auto collisionsPred = [](Container_1 a, const auto &b) {

        for( auto& c : a) {
            // Check if any of the collisions objects are in the collection already
            if( c.obj1 == b.obj1 )
                return true;
            if( c.obj1 == b.obj2 )
                return true;
            if( c.obj2 == b.obj1 )
                return true;
            if( c.obj2 == b.obj2 )
                return true;

            return false;// add a new object to container of collisions
        }
    };

    auto singularitiesPred = [](Container_2 a, const auto &b) {

        for( auto& c : a) {

            // Check if any of the state objects are in the collection already
            if( c.obj == b.obj ) return true; //for singularities

            return false;
        }
    };

    auto start = std::end(container2);

    // Check if any object that is in CONTAINER1 is also in CONTAINER2 and compares time values to remove one
    for( auto first_iter = std::end(container1) - 1; first_iter != std::begin(container1) -1; --first_iter) {
        for( auto second_iter = start - 1; second_iter != std::begin(container2) -1; --second_iter ) {


            if(( objPred( *first_iter, *second_iter ) )) { //check if the same object

                if( timePred( *first_iter, *second_iter )) { //compare the time

                    // Keep collision
                    if( collisionsPred(realColContainer, *first_iter) == false) {//if collision first

                        realColContainer.push_back(*first_iter);
                    }
                }
                else {

                    // Keep state
                    if( singularitiesPred(realSingContainer, *second_iter) == false ) {//if singularity is first

                        realSingContainer.push_back(*second_iter);
                    }
                }
            }
        }
    }

    _collisions = realColContainer;
    _singularities = realSingContainer;
}
#define end1 }


#define sortAndMakeUniqueStatesMethod {
// Sort and make Unique for states
template <class Container>
void MyController::sortAndMakeUniqueStates(Container& c) {

    // Sort
    std::sort( std::begin(c), std::end(c), [] (const auto& a, const auto&b) {
        return a.time < b.time;
    });

    // Make unique
    auto pred =  [](const auto &a, const auto &b) {

        if( a.obj == b.obj ) return true;

        return false;
    };

    typename Container::iterator NE = std::end(c);
    for( auto first_iter = std::begin(c); first_iter != NE; ++first_iter) {

        for( auto r_iterator = NE - 1; r_iterator != first_iter; --r_iterator) {

            if( (pred(*first_iter, *r_iterator))) {
                std::swap( *r_iterator, *(NE-1) );
                NE--;

            }
        }
    }

    c.erase( NE, std::end( c ) );
}
#define end1 }

#define sortAndMakeUniqueCollisionsMethod {
template <class Container_T >
void sortAndMakeUnique( Container_T& container) {

    // Sort
    std::sort( std::begin(container), std::end(container), [](const auto& a, const auto& b) {
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

}



#define detectCollisionsMethod {
// Sphere and Static Plane
template <typename T_s>
inline
void MyController::sphereCollision(T_s* sphere, seconds_type dt)
{
    bool moveObj = false;
    //check if any of spheres are moving
    for( auto& sphere : _dynamic_spheres ) {
        if( sphere->_state != States::Still ) {
            moveObj = true;
            break;
        }
    }

    if( moveObj ) {

        //start
        for( auto iter1 = std::begin(_dynamic_spheres); iter1 != std::end(_dynamic_spheres); ++iter1) {


            auto col = detectCollision(*sphere,**iter1,dt);

            if(col.CollisionState::flag == CollisionStateFlag::Collision){
                auto& sphere1 = sphere;
                auto& sphere2 = *iter1;

                auto new_dt = std::max(sphere1->curr_t_in_dt,sphere2->curr_t_in_dt);

                if (col.time > seconds_type(new_dt) and col.time < dt){

                    auto col_obj = CollisionObject(sphere1,sphere2,col.time);
                    _collisions.push_back((col_obj));

                }
            }
        }


        for( auto iter1 = std::begin(_static_planes); iter1 != std::end(_static_planes); ++iter1) {


            if (sphere->_state != States::Still){

                auto col = detectCollision(*sphere,**iter1,dt);

                if(col.CollisionState::flag == CollisionStateFlag::Collision){

                    auto new_dt = sphere->curr_t_in_dt;

                    if (col.time > seconds_type(new_dt) and col.time < dt){

                        auto col_obj = CollisionObject(sphere,*iter1,col.time);
                        _collisions.push_back((col_obj));

                    }
                }

            }
        }
    }

}
#define end1 }

} // END namespace collision



#endif





