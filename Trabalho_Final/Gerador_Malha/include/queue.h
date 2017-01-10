#ifndef QUEUE_H
#define QUEUE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

// =============================================================================================================
// Estrutura de um Node da fila
struct QNode
{
	QNode *next;				// Ponteiro para o proximo elemento da fila
	QNode *prev;				// Ponteiro para o elemento anterior da fila
	int node;					// Indice do nodo em crescimento da arvore
}typedef QNode;

QNode* newQNode (QNode *next, QNode *prev, int node);
// =============================================================================================================
// =============================================================================================================
// Estrutura da fila
struct Queue
{
	QNode *head;				// Ponteiro para a cabeca da fila
	QNode *last;				// Ponteiro para o ultimo da fila
	int in_the_queue;			// Contador de nos em crescimento na fila
}typedef Queue; 

Queue* newQueue ();
void Enqueue (Queue **q, int p);
bool isEmpty (Queue *q);
int Dequeue (Queue **q);
void printQueue (Queue *q);
// =============================================================================================================
// =============================================================================================================

#endif