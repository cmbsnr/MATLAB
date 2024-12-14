%% Neutral DDE with Two Delays
% Solve the following neutral DDE, presented by Paul, for $0 \le t \le
% \pi$.
%
% $$y'(t) = 1 + y(t) - 2y(t/2)^2 - y'(t-\pi)$$
%
% The solution history is $y(t) = \cos(t)$ for $t \le 0$.

%%
% Create a new program file in the editor. This file will contain a main
% function and four local functions.

%%
% Define the first-order DDE as a local function named |ddefun|.
%
% <include>ddefun.m</include>
%

%%
% Define the solution delay as a local function named |dely|.
%
% <include>dely.m</include>
%

%%
% Define the derivative delay as a local function named |delyp|.
%
% <include>delyp.m</include>
%

%%
% Define the solution history as a local function named |history|.
%
% <include>history.m</include>
%

%%
% Define the interval of integration and solve the DDE using |ddensd|. Add
% this code to the main function.
tspan = [0 pi];
sol = ddensd(@ddefun,@dely,@delyp,@history,tspan);

%%
% Evaluate the solution at 100 equally spaced points between $0$ and $\pi$.
% Add this code to the main function.
tn = linspace(0,pi);
yn = deval(sol,tn);

%%
% Plot the results. Add this code to the main function.
plot(tn,yn);
xlim([0 pi]);
ylim([-1.2 1.2]);
xlabel('time t');
ylabel('solution y');

%%
% Run your entire program to calculate the solution and display the plot.
% The file |ddex4.m| contains the complete code for this example. To see
% the code in an editor, type |edit ddex4| at the command line.


%% 
% Copyright 2012 The MathWorks, Inc.