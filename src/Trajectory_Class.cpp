#include "Trajectory_Class.hpp"

#include <iostream>

#include "Physical_Parameters.hpp"
#include "General_Utilities.hpp"

//Event Class
	//Constructors
		Event::Event(){
			time=0.0;
			position=Eigen::Vector3d (0.0,0.0,0.0);
			velocity=Eigen::Vector3d (0.0,0.0,0.0);
		}
		Event::Event(double t, Eigen::Vector3d& x, Eigen::Vector3d& v){
			time=t;
			position=x;
			velocity=v;
		}
		Event::Event(double t, double x, double y, double z, double vx, double vy, double vz){
			time=t;
			position=Eigen::Vector3d (x,y,z);
			velocity=Eigen::Vector3d (vx,vy,vz);
		}
	//Set Values
		void Event::SetTime(double t){
			time=t;
		}
		void Event::SetPosition(Eigen::Vector3d& newposition){
			position=newposition;
		}
		void Event::SetVelocity(Eigen::Vector3d& newvelocity){
			velocity=newvelocity;
		}
		void Event::SetPosition(double x,double y,double z){
			position=Eigen::Vector3d (x,y,z);
		}
		void Event::SetVelocity(double vx,double vy,double vz){
			velocity=Eigen::Vector3d (vx,vy,vz);
		}
	//Change values
		void Event::IncreaseTime(double dt)
		{
			time+=dt;
		}
		void Event::IncreasePosition(Eigen::Vector3d& dx)
		{
			position+=dx;
		}
	//Return Values
		double Event::Time()
		{
			return time;
		}
		Eigen::Vector3d Event::Position(){
			return position;
		}
		Eigen::Vector3d Event::Velocity()
		{
			return velocity;
		}
		//Returns the event in sec, km, km/sec.
		Event Event::kmsec()
		{
			double t=time/sec;
			Eigen::Vector3d x=1/km*position;
			Eigen::Vector3d v=sec/km*velocity;
			return Event(t,x,v);
		}
	//Vector Norms
		double Event::NormPosition(){
			return position.norm();
		}
		double Event::NormVelocity(){
			return velocity.norm();
		}
	//Overload <<
		std::ostream& operator<<(std::ostream &output,const Event &event){
			return output 	<<"{"
							<<event.time
							<<","
							<<event.position.format(VectorFormat)
							<<","
							<<event.velocity.format(VectorFormat)
							<<"}";
		}



//Trajectory Class
	//Constructors:
		Trajectory::Trajectory(){
			weight=1.0;
		}
		Trajectory::Trajectory(std::vector<Event> input){
			weight=1.0;
			events=input;
		}
		Trajectory::Trajectory(std::vector<Event> input,double w)
		{
			weight=w;
			events=input;
		}
	//Extract Info from a single Track
		int Trajectory::Length()
		{
			return events.size();
		}

		double Trajectory::tStart(){
			return events.front().Time();
		}
		double Trajectory::tEnd(){
			return events.back().Time();
		}
		double Trajectory::tEarthEntry(){
			return events[1].Time();
		}
		Event Trajectory::TrackInterpolation(double t){
			if(t<tStart()||t>tEnd())
			{
				std::cout <<"Error: TrackInterpolation was given an argument outside the Trajectorys range."<<std::endl;
				printf("Range is [%f,%f], but t = %f\n",tStart()/sec,tEnd()/sec,t);
				return Event();
			}
			else
			{
				for(unsigned int i=1;i<events.size();i++)
				{
					if(t<=events[i].Time())
					{
						Eigen::Vector3d x=events[i-1].Position()+(t-events[i-1].Time())*events[i-1].Velocity();
						Eigen::Vector3d v=events[i-1].Velocity();
						return Event(t,x,v);
					}
				}
			}
			cout <<"Error in TrackInterpolation(): This shouldn't happen..."<<endl;
			return Event();	
		}
		Event Trajectory::TrajectoryStart()
		{
			return events.front();
		}
		Event Trajectory::TrajectoryEnd()
		{
			return events.back();
		}

		void Trajectory::New_Event(Event event,double weight)
		{
			Update_Weight(weight);
			events.push_back(event);
		}
		//We differentiate three types of Trajectory.
		// 1 - particle enters the earth and leaves
		// 2 - particle enters the earth and does not leave (typically, the velocity cutoff is reached)
		// 3 - particle reaches the detector depth  
		int Trajectory::Trajectory_Type()
		{
			int Wcase;
			if(events.back().Position()[2]>0) Wcase=1;
			else 
			{
				
				//particle reaches detector depth
				if( abs( events.back().Position()[2] + Detector_Depth ) < cm  )
				{
					Wcase =3;
				}
				//particle gets stuck case 2
				else
				{
					Wcase =2;
				}
			}
			return Wcase;
		}

		double Trajectory::Final_Speed()
		{
			return TrajectoryEnd().Velocity().norm();
		}

		int Trajectory::NoOfScatterings()
		{
			int nScatterings;
			int Wcase =  Trajectory_Type();
			if (Wcase==1) nScatterings=events.size()-4;
			else if (Wcase==2) nScatterings=events.size()-2;
			else if (Wcase==3)nScatterings=events.size()-3;
			else nScatterings=0;
			return nScatterings;
		}
		double Trajectory::TrajectoryLength()
		{
			double length=0.0;
			for(unsigned int i=1;i<events.size();i++)
			{
				length+= (events[i].Position()-events[i-1].Position()).norm();
			}
			return length;
		}

		double Trajectory::MaximumDepth()
		{
			double maxDepth=0.0;
			for(unsigned int i=0;i<events.size();i++)
			{
				if(events[i].Position()[2]<maxDepth) maxDepth=events[i].Position()[2];
			}
			return (-1.0)*maxDepth;
		}

		double Trajectory::Maximum_Horizontal_Distance()
		{
			double lMax =0.0;
			for(unsigned int i = 0; i<events.size();i++)
			{
				double l = sqrt( pow(events[i].Position()[0] - events[0].Position()[0],2.0) + pow(events[i].Position()[0] - events[0].Position()[0],2.0));
				if(l>lMax) lMax=l;
			}
			return lMax;
		}


		double Trajectory::Get_Weight()
		{
			return weight;
		}

		void Trajectory::Update_Weight(double w)
		{
			weight*=w;
		}

		std::vector<Event> Trajectory::DepthCrossing(double d)
		{
			vector<Event> EventList;

			return EventList;
		}

	//Overload <<
		std::ostream& operator<<(std::ostream &output,const Trajectory &trajectory){
			output <<"{";
			for(unsigned int i=0;i<trajectory.events.size();i++)
			{
				output << trajectory.events[i];
				if(i<(trajectory.events.size()-1))
				{
					output <<",";
				}
			}
			output <<"}";
			return output;
		}
		
