#include <iostream>

#include "add.h"
#include "subtract.h"

int main()
{
    double a = 10;
    double b = 3;

    std::cout << "Add: " << a << " + " << b
              << " = " << add(a, b) << std::endl;

    std::cout << "Subtract: " << a << " - " << b
              << " = " << subtract(a, b) << std::endl;

    system("PAUSE");
    
    return 0;
}
