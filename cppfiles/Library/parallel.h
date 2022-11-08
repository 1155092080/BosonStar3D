#include <thread>
using namespace std;

template <typename Fun, typename... Args>
void parallel_for(int start, int end, Fun fun, Args... args){
	
	int noOfElement = end - start;
	int noOfThread = thread::hardware_concurrency();
	
	if(noOfElement > noOfThread){
		thread threads[noOfThread];
		for(int i=0;i<noOfThread;i++){
			threads[i] = thread(fun, start+i, args...);
		}
		
		for(int i=0;i<noOfThread;i++){
			threads[i].join();
		}
		
		parallel_for(start+noOfThread, end, fun, args...); 
		
	}else{
		thread threads[noOfElement];
		for(int i=0;i<noOfElement;i++){
			threads[i] = thread(fun, start+i, args...);
		}
		
		for(int i=0;i<noOfElement;i++){
			threads[i].join();
		}
	}
	
	return;
	
}

