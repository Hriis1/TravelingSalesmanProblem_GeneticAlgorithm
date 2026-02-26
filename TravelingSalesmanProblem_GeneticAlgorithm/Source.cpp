#include <iostream>
#include <vector>

void outputPath(const std::vector<int>& path)
{
	const int sz = path.size();

	if (sz == 0)
		return;

	//Output path
	for (size_t i = 0; i < sz; i++)
	{
		std::cout << path[i] << " -> ";
	}

	//Return to start city
	if (sz > 1) 
	{
		std::cout << path[0];
	}

	//New line
	std::cout << std::endl;
}

int main() 
{
	std::cin.get();
}