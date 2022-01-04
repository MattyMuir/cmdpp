#pragma once
#include <stdio.h>
#include <iostream>
#include <vector>
#include <memory>

#include <mpir.h>

class mpzxx
{
public:
	~mpzxx()
	{
		mpz_clear(data);
	}
	mpz_t data;
};

class mpzArray
{
public:
	// Functions
	mpzArray();
	mpzArray(int count);
	~mpzArray();

	void Reserve(int count);
	void Clear();
	void TempClear();
	mpz_t& operator[](int i);
	int Size();
	void Back();
	void Back(mpz_t element);
	void BackUI(uint64_t element);
	void DeepCopy(mpzArray& other);
	void Resize(int newSize);

private:
	std::vector<std::shared_ptr<mpzxx>> mData;
	int size;
};