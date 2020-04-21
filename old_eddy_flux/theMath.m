clear;close all;clc

syms Top a To Tap Ta U Up CD

expand(CD*(1+a*Top)*(To+Top-Ta-Tap)*(U+Up))

CD*To*U - CD*Tap*U - CD*Ta*U + CD*Top*U -...
CD*Ta*Up - CD*Tap*Up + CD*To*Up +...
CD*Top*Up + CD*Top^2*U*a + CD*Top^2*Up*a -...
CD*Ta*Top*U*a - CD*Tap*Top*U*a + CD*To*Top*U*a -...
CD*Ta*Top*Up*a - CD*Tap*Top*Up*a + CD*To*Top*Up*a

% spatial average

CD*To*U - CD*Ta*U  - CD*Tap*Up  + CD*Top*Up +
CD*Top^2*U*a + CD*Top^2*Up*a  - CD*Tap*Top*U*a 
- CD*Ta*Top*Up*a - CD*Tap*Top*Up*a + CD*To*Top*Up*a






