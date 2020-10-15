#include <iostream>
#include <vector>

using namespace std;

int main(){
	const int big_size = 10000;
	vector<double> v( big_size );
	cout << "Before clearing, the capacity of the vector is "
		  << v.capacity() << " and its size is " << v.size();
	v.clear();
	cout << "\nAfter clearing, the capacity of the vector is "
		  << v.capacity() << " and its size is " << v.size();
	vector<double>().swap( v );

	cout << "\nAfter swapping, the capacity of the vector is "
		  << v.capacity() << " and its size is " << v.size();

	vector<double> v1( big_size );
	v1 = vector<double>();
	cout << "\n After vector=<double>();, the capacity of the vector is "
		  << v1.capacity() << " and its size is " << v1.size() << endl;
	return 0;
}
