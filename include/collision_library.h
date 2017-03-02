#ifndef COLLISION_LIBRARY_H
#define COLLISION_LIBRARY_H

// collision library interface
#include <collision_interface.h>

#include <typeinfo>
#include <unordered_map>


namespace collision
{


//States
enum class States {
    Free,
    Rolling,
    Still,
};


//struct StateChangeObj;
struct StateChangeObj {
    DynamicPSphere*                             obj;               // Object whos state will change
    std::unordered_set<StaticPPlane*>           attachedPlanes;     // Object that obj1 will change state ACCORDING to
    seconds_type                                time;               // Time of singularity
    States                                      state;              // State that obj1 will change to

    StateChangeObj
    (DynamicPSphere* o1, std::unordered_set<StaticPPlane*> planes, seconds_type t, States s) :
        obj{o1}, attachedPlanes{planes}, time{t} ,state{s} {}
};


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
    //    std::vector<StaticPPlane*> getAttachedObjects(DynamicPSphere* sphere);


    void setAttachedObjects( std::unordered_set<StaticPPlane*> plane, DynamicPSphere* sphere );
    void detachObjects(DynamicPSphere* sphere);


    void
    detectStateChanges(double dt);//fill the container of singularities

    StateChangeObj
    detectStateChange(DynamicPSphere* sphere, double dt); //check the state

    void
    handleStates (StateChangeObj& state, double dt);//update the state

    GMlib::Point<float,2>
    closestPoint(DynamicPSphere& S, StaticPPlane& P);//get closest point

    void handleCollision ( collision::CollisionObject& col, double dt);//checks what objects are collided

    Environment                                 _stillEnvironment;
    DefaultEnvironment                          _env;

    std::unordered_map<DynamicPSphere*, std::unordered_set<StaticPPlane*>>    _attachedPlanes;


private:

    template <typename T_s, typename T_o>
    void sphereStaticCollision(T_s* sphere, T_o* object, seconds_type dt);//collison detection

    template <typename T_s, typename T_o>
    void sphereDynamicCollision(T_s* sphere, T_o* object, seconds_type dt);

    template <class Container>
    void sortAndMakeUniqueStates(Container& c);//sort and make states

    template <typename Container_1, typename Container_2>
    void crossUnique( Container_1 container1, Container_2 container2);// take two sorted and unique containers and crosscheck them,
    //update containers

protected:
    void localSimulate (double dt) override;

    std::vector<DynamicPSphere*>                _dynamic_spheres;
    std::vector<StaticPSphere*>                 _static_spheres;
    std::vector<StaticPPlane*>                  _static_planes;
    std::vector<StaticPCylinder*>               _static_cylinders;
    std::vector<StaticPBezierSurf*>             _static_bezier_surf;

    std::vector<collision::CollisionObject>     _collisions;
    std::vector<StateChangeObj>                 _singularities;





};

template <>
class DynamicPhysObject<GMlib::PSphere<float>> : public DynamicPhysObject_Base<GMlib::PSphere<float>> {
public:
    using DynamicPhysObject_Base<GMlib::PSphere<float>>::DynamicPhysObject_Base;


    MyController*   _sphereController;
    States          _state = States::Free;     // Which state is the sphere in

    bool            checker = false;

    void moveUp();
    void moveDown();
    void moveLeft();
    void moveRight();

    float                       u;
    float                       v;


    GMlib::Vector<float,3> getSurfNormal();
    void setUV(StaticPPlane* plane);

    void computeStep(double dt);

    void updateX(double x);
    double getX();


    void setVelocity(const GMlib::Vector<double,3> vel);
    GMlib::Vector<double,3> getVelocity();

    double getMass();
    GMlib::Vector<float,3> getDs();


    GMlib::Vector<double,3>
    adjustedTrajectory (seconds_type dt);


    void            simulateToTInDt( seconds_type t ) override;

    GMlib::Vector<double, 3>
    computeTrajectory (seconds_type dt) const override; // [m]

    GMlib::Vector<double, 3> externalForces () const override; // [m / s^2]


};


template <class PSurf_T, typename... Arguments>
std::unique_ptr<DynamicPhysObject<PSurf_T>> unittestDynamicPhysObjectFactory(Arguments... parameters) {

    return std::make_unique<DynamicPhysObject<PSurf_T>> (parameters...);

}

template <class PSurf_T, typename... Arguments>
std::unique_ptr<StaticPhysObject<PSurf_T>> unittestStaticPhysObjectFactory(Arguments... parameters) {

    return std::make_unique<StaticPhysObject<PSurf_T>> (parameters...);
}

// Make collisions and states crosswise unique
template <typename Container_1, typename Container_2>
void MyController::crossUnique(Container_1 container1, Container_2 container2) {

    auto objPred =  [](const auto &a, const auto &b) {

        if( a.obj1 == b.obj or a.obj2 == b.obj ) return true; //check if equal

        return false;
    };

    auto timePred =  [](const auto &a, const auto &b) {

        if( a.t_in_dt < b.time ) return true; //compare objects time, if collision is before the state or not

        return false;//if state is first
    };

    auto inCollisionsAlreadyPred = [](Container_1 a, const auto &b) {

        for( auto& c : a) {

            // Check if any of the collisions objects are in the collection already
            if( c.obj1 == b.obj1 ) return true;
            if( c.obj1 == b.obj2 ) return true;
            if( c.obj2 == b.obj1 ) return true;
            if( c.obj2 == b.obj2 ) return true;

            return false;// add a new object to container of collisions
        }
    };

    auto inSingularitiesAlreadyPred = [](Container_2 a, const auto &b) {

        for( auto& c : a) {

            // Check if any of the state objects are in the collection already
            if( c.obj == b.obj ) return true; //for singularities

            return false;
        }
    };

    std::vector<collision::CollisionObject>     _realCollisions; //new containers
    std::vector<StateChangeObj>                 _realSingularities;

    auto start = std::end(container2);

    // Check if any object that is in CONTAINER1 is also in CONTAINER2 and compares time values to remove one
    for( auto first_iter = std::end(container1) - 1; first_iter != std::begin(container1) -1; --first_iter) {
        for( auto second_iter = start - 1; second_iter != std::begin(container2) -1; --second_iter ) {


            if(( objPred( *first_iter, *second_iter ) )) { //check if the same object

                if( timePred( *first_iter, *second_iter )) { //compare the time

                    // Keep collision
                    if( inCollisionsAlreadyPred(_realCollisions, *first_iter) == false) {//if collision first

                        _realCollisions.push_back(*first_iter);
                    }
                }
                else {

                    // Keep state
                    if( inSingularitiesAlreadyPred(_realSingularities, *second_iter) == false ) {//if singularity is first

                        _realSingularities.push_back(*second_iter);
                    }
                }
            }
        }
    }

    _collisions = _realCollisions;
    _singularities = _realSingularities;
}

// Sort and make Unique for states
template <class Container>
void MyController::sortAndMakeUniqueStates(Container& c) {

    // Sorting
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

} // EOF


template <typename T_s, typename T_o>
inline
void MyController::sphereDynamicCollision(T_s* sphere, T_o* object, seconds_type dt) {

    bool movingObject = false;

    for( auto& s : _dynamic_spheres ) {
        if( s->_state != States::Still ) {
            movingObject = true;
            break;
        }
    }

    if( movingObject ) {

        // Sphere vs. Dynamic Sphere
        for( auto iter = std::begin(_dynamic_spheres); iter != std::end(_dynamic_spheres); ++iter) {

            if( (*iter == sphere) or (*iter == object) or (sphere->_state == States::Still)) {//if sphere is not the same or still
                break;
            }
            else {

                // Check only for collisions for the first dynamic sphere
                auto col = collision::detectCollision(*sphere, **iter, dt);

                // The check for if the second object is a dynamic sphere is done in the main algorithm
                // and the method is called again only with the 1'st and 2'nd sphere swapped

                if( col.CollisionState::flag == CollisionStateFlag::Collision ) {

                    auto& first_sphere = sphere;
                    auto& second_sphere = *iter;

                    auto new_t = std::max(first_sphere->curr_t_in_dt, second_sphere->curr_t_in_dt);

                    if( col.time > seconds_type(new_t) and col.time < dt) {
                        auto col_obj = CollisionObject(first_sphere, second_sphere, col.time);
                        _collisions.push_back(col_obj);


                    }
                }
            }
        }
    }

}

template <typename T_s, typename T_o>
inline
void MyController::sphereStaticCollision(T_s* sphere, T_o* object, seconds_type dt)
{

    bool movingObject = false;

    for( auto& s : _dynamic_spheres ) {
        if( s->_state != States::Still ) {
            movingObject = true;
            break;
        }
    }

    if( movingObject ) {


        // The game will mainly have two types of static objects that can be collided with, planes and bezier surfaces
        // Checks to see which one the sphere has collided with
        auto plane = dynamic_cast<const StaticPPlane*>(object);
        //auto bezier_s = dynamic_cast<const StaticPBezierSurf*>(object);


        // Sphere vs. Static Plane
        if( plane ) {

            auto attachedPlanes = getAttachedObjects(sphere);

            for( auto iter1 = std::begin(_static_planes); iter1 != std::end(_static_planes); ++iter1) {


                // If sphere is attached to the plane, it should NOT check for collisions with it
                if( !attachedPlanes.empty() ) {
                    for( auto iter2 = attachedPlanes.begin(); iter2 != attachedPlanes.end(); ++iter2) {

                        if( *iter1 == *iter2) break;
                        else {

                            // Else look for collision
                            auto col = collision::detectCollision(*sphere, **iter1, dt);

                            if( col.CollisionState::flag == CollisionStateFlag::Collision ) {

                                auto& static_object = *iter1;

                                auto new_t = sphere->curr_t_in_dt;

                                if( col.time > seconds_type(new_t) and col.time < dt)
                                {
                                    auto col_obj = CollisionObject(sphere, static_object, col.time);
                                    _collisions.push_back(col_obj);

                                }
                            }
                        }
                    }
                }
                else if( sphere->_state != States::Still) {

                    auto col = collision::detectCollision(*sphere, **iter1, dt);

                    if( col.CollisionState::flag == CollisionStateFlag::Collision ) {

                        auto& static_object = *iter1;

                        auto new_t = sphere->curr_t_in_dt;

                        if( col.time > seconds_type(new_t) and col.time < dt)
                        {
                            auto col_obj = CollisionObject(sphere, static_object, col.time);
                            _collisions.push_back(col_obj);

                        }
                    }
                }
            }
        }

    }


}


} // END namespace collision



#endif





