#include <iostream>
using namespace std;
int main()
{
    

    int x = 2;
    int* ptr = &x;

    cout<<x<<endl;
    //prints 2

    cout<<ptr<<endl;
    //prints address of x

    cout<<*ptr<<endl;
    //prints 2
}