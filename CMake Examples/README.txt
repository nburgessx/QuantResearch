How to generate the Windows Solution File?
------------------------------------------

1. Open Command Line (with cmake)
From Visual Studio open a Visual Studio Command Prompt

2. Goto Build Folder
Navigate to the build folder of the example project

3. Generate the Windows Solution File using CMake
On command line type:
cmake -G "Visual Studio 17 2022" ..

Generally we type:
cmake -G "Compiler Name" %(path to root cmake folder)

The switch -G is case sensitive and refers to the generator
Here we use '..' as the root folder is one level above the build folder

To build the code on linux change the compiler name e.g.
cmake -G "Unix Makefiles" ..
cmake -G Ninja ...
Note: We use quotes if the compiler name contains a space.
