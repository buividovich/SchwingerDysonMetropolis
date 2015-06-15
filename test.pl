#!/usr/bin/perl

@A = (1, 2);

$x1 = pop(@A); $nA = scalar(@A);
print("x1 = $x1, A after pop: @A, scalar(A) = $nA\n");

$x2 = pop(@A); $nA = scalar(@A);
print("x2 = $x2, A after pop: @A, scalar(A) = $nA\n");

