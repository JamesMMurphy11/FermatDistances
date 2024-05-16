function uTest_Cell2Vec(doSpeed)
% Automatic test: Cell2Vec
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_Cell2Vec(doSpeed)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2009-2015 matlab.2010(a)n(MINUS)simon.de

% $JRev: R-u V:020 Sum:mX6sLlK2RFdj Date:27-Aug-2012 00:38:38 $
% $License: BSD (see Docs\BSD_License.txt) $
% $File: Tools\UnitTests_\uTest_Cell2Vec.m $
% History:

% Initialize: ==================================================================
ErrID = ['JSimon:', mfilename];
LF = char(10);

% Default for the input:
if nargin == 0
   doSpeed = true;
end

% Times for testing:
if doSpeed
   randDelay = 5;    % [sec], time for random tests
   SpeedTime = 1.0;  % [sec], time for speed tests
else
   randDelay = 0.5;
   SpeedTime = 0.1;
end
TestTime = 0.2;  % Get number of loops

% Hello:
whichFunc = which('Cell2Vec');
disp(['==== Test Cell2Vec:  ', datestr(now, 0), LF, ...
      'Version: ', whichFunc, LF]);

% Determine if the MEX or M-version is called:
[dummy1, dummy2, fileExt] = fileparts(whichFunc);  %#ok<ASGLU>
useMex = strcmpi(strrep(fileExt, '.', ''), mexext);

% Start tests: -----------------------------------------------------------------
% Standard tests - empty input, cell without populated elements, small known
% answer test:
S = Cell2Vec({});
if isempty(S) && isa(S, 'double')
   disp('  ok: empty cell');
else
   error(ErrID, 'Failed for empty cell.');
end

S = Cell2Vec(cell(1, 1));
if isempty(S) && isa(S, 'double')
   disp('  ok: {NULL}');
else
   error(ErrID, 'Failed for {NULL}.');
end

S = Cell2Vec(cell(1, 2));
if isempty(S) && isa(S, 'double')
   disp('  ok: {NULL, NULL}');
else
   error(ErrID, 'Failed for {NULL, NULL}.');
end

S = Cell2Vec({''});
if isempty(S) && ischar(S)
   disp('  ok: {''''}');
else
   error(ErrID, 'Failed for {''''}.');
end

S = Cell2Vec({'', ''});
if isempty(S) && ischar(S)
   disp('  ok: {'''', ''''}');
else
   error(ErrID, 'Failed for {'''', ''''}.');
end

S = Cell2Vec({'S'});
if isequal(S, 'S')
   disp('  ok: {''S''}');
else
   error(ErrID, 'Failed for {''S''}.');
end

S = Cell2Vec({'S1', 'S2'});
if isequal(S, 'S1S2')
   disp('  ok: {''S1'', ''S2''}');
else
   error(ErrID, 'Failed for {''S1'', ''S2''}.');
end

S = Cell2Vec({'S'; 'T'; 'UV'});
if isequal(S, 'STUV')
   disp('  ok: {''S''; ''T''; ''UV''}');
else
   error(ErrID, 'Failed for {''S''; ''T''; ''UV''}.');
end

S = Cell2Vec({'AAA', 'BB', [], 'C'});
if isequal(S, 'AAABBC')
   disp('  ok: {''AAA''; ''BB''; [], ''C''}');
else
   error(ErrID, 'Failed for {''AAA''; ''BB''; [], ''C''}.');
end

x = rand(10);
y = rand(11);
S = Cell2Vec({x, y});
if isequal(S, transpose([x(:); y(:)]))
   disp('  ok: {rand(10), rand(11)}');
else
   error(ErrID, 'Failed for {rand(10), rand(11)}.');
end

x = single(x);
y = single(y);
S = Cell2Vec({x, y});
if isequal(S, transpose([x(:); y(:)]))
   disp('  ok: {single(rand(10)), single(rand(11))}');
else
   error(ErrID, 'Failed for {single(rand(10)), single(rand(11))}.');
end

S = Cell2Vec({true, false});
if isequal(S, [true, false])
   disp('  ok: {TRUE, FALSE}');
else
   error(ErrID, 'Failed for type LOGICAL.');
end

if sscanf(version, '%d', 1) >= 7.0
   TypeList = {'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', ...
         'int64', 'uint64'};
else  % No RESHAPE or TRANSPOSE for UINT64 in Matlab 6.5...
   TypeList = {'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32'};
end

for iType = 1:length(TypeList)
   aType = TypeList{iType};
   x = feval(aType, round(rand(3,5) * 127));
   y = feval(aType, round(rand(3,5) * 127));
   S = Cell2Vec({x, y});
   if isequal(S, reshape([x(:); y(:)], 1, []))
      disp(['  ok: {', aType, ', ', aType, '}']);
   else
      error(ErrID, ['Failed for type ', upper(aType)]);
   end
end

% Random tests: ----------------------------------------------------------------
fprintf('\n== Random tests (%g sec):\n', randDelay);
lenDataC = 100;
DataC    = cell(1, lenDataC);
for k = 1:lenDataC
   DataC{k} = sprintf('%d', fix(1000 .^ rand));
end

iniTime = cputime;
nTest   = 0;
while cputime - iniTime < randDelay
   for N = 0:99
      C = DataC(fix(rand(1, N) * lenDataC) + 1);
      if ~isequal(horzcat(C{:}), Cell2Vec(C))
         error(ErrID, 'Failed for random test.');
      end
   end
   nTest = nTest + 100;
end
fprintf('  ok: %d random tests passed.\n', nTest);

% Invalid input: ---------------------------------------------------------------
fprintf('\n== Check rejection of bad input:\n');
tooLazy = false;

try
   dummy   = Cell2Vec([]);
   tooLazy = true;
catch
   disp(['  ok: [] rejected: ', LF, '      ', strrep(lasterr, LF, '; ')]);
end
if tooLazy
   error(ErrID, '[] not rejected.');
end

if useMex  % Less checks of inputs in M-version:
   try
      dummy   = Cell2Vec({1, 'a'});  %#ok<*NASGU>
      tooLazy = true;
   catch
      disp(['  ok: {1, ''a''} rejected: ', LF, ...
            '      ', strrep(lasterr, LF, '; ')]);
   end
   if tooLazy
      error(ErrID, '{1, ''a''} not rejected.');
   end
   
   try
      dummy   = Cell2Vec({'1', 2});
      tooLazy = true;
   catch
      disp(['  ok: {''1'', 2} rejected: ', LF, '      ', ...
            strrep(lasterr, LF, '; ')]);
   end
   if tooLazy
      error(ErrID, '{''1'', 2} not rejected.');
   end
   
   try
      dummy   = Cell2Vec({'1'}, []);
      tooLazy = true;
   catch
      disp(['  ok: ({''1''}, []) rejected: ', LF, '      ', ...
            strrep(lasterr, LF, '; ')]);
   end
   if tooLazy
      error(ErrID, '({''1''}, []) not rejected.');
   end
else
   disp('  ?:  Reduced input checks for M-version');
end

% Speed: -----------------------------------------------------------------------
if doSpeed
   fprintf('\n== Speed test (test time: %g sec):\n', SpeedTime);
else
   fprintf('\n== Speed test (test time: %g sec - may be inaccurate):\n', ...
      SpeedTime);
end
drawnow;  % Allow update of external events

% TypeList = {'int8', 'int16', 'single', 'double'};
TypeList = {'double'};
CellSize = {[1,5], [1,10], [1,100], [10,1000]};
DataSize = {[1,1], [1,10], [1,100], [1,1000]};
fprintf('  Cell size:           ');
for i = 1:length(CellSize)
   fprintf('%9s', sprintf('%dx%d', CellSize{i}));
end
fprintf('\n');
rejectLargeData = false;
for iType = 1:length(TypeList)
   aType = TypeList{iType};
  
   for iDataSize = 1:length(DataSize)
      aDataSize = DataSize{iDataSize};
      fprintf('  Data: [%d x %d] %s\n', aDataSize, aType);
      
      fprintf('    CAT(2, C{:})       ');
      for iCellSize = 1:length(CellSize)
         aCellSize = CellSize{iCellSize};
         if rejectLargeData && ...
               prod(aDataSize) >= 10000 && prod(aCellSize) >= 10000
            PrintLoop([]);
            continue;
         end
         
         CStr = cell(aCellSize);
         for iC = 1:prod(aCellSize)
            CStr{iC} = feval(aType, rand(aDataSize) * 255);
         end
         
         % Get number of loops (the WHILE CPUTIME loop has more overhead):
         iTime = cputime;
         N = 0;
         while cputime - iTime < TestTime
            S = cat(2, CStr{:});
            clear('S');
            N = N + 1;
         end
         N = max(2, round(N / (cputime - iTime) * SpeedTime));
         
         % The actual test:
         tic;
         for iN = 1:N
            S = cat(2, CStr{:});
            clear('S');
         end
         NPerSec = N / (toc + eps);  % Loops per second
         PrintLoop(NPerSec);
         drawnow;  % Allow update of external events
      end
      fprintf('  loops/sec\n');
      
      % ********
      fprintf('    [C{:}]  (HORZCAT!) ');
      for iCellSize = 1:length(CellSize)
         aCellSize = CellSize{iCellSize};
         if rejectLargeData && ...
               prod(aDataSize) >= 10000 && prod(aCellSize) >= 10000
            PrintLoop([]);
            continue;
         end
         
         CStr = cell(aCellSize);
         for iC = 1:prod(aCellSize)
            CStr{iC} = feval(aType, rand(aDataSize) * 255);
         end
         
         % Get number of loops (the WHILE CPUTIME loop has more overhead):
         iTime = cputime;
         N = 0;
         while cputime - iTime < TestTime
            S = [CStr{:}];
            clear('S');
            N = N + 1;
         end
         N = max(2, round(N / (cputime - iTime) * SpeedTime));
         
         % The actual test:
         tic;
         for iN = 1:N
            S = [CStr{:}];
            clear('S');
         end
         NPerSec = N / (toc + eps);  % Loops per second
         PrintLoop(NPerSec);
         drawnow;  % Allow update of external events
      end
      fprintf('  loops/sec\n');
      
      % ********
      fprintf('    CELL2MAT(C)        ');
      for iCellSize = 1:length(CellSize)
         aCellSize = CellSize{iCellSize};
         if rejectLargeData && ...
               prod(aDataSize) >= 10000 && prod(aCellSize) >= 10000
            PrintLoop([]);
            continue;
         end
         
         CStr = cell(aCellSize);
         for iC = 1:prod(aCellSize)
            CStr{iC} = feval(aType, rand(aDataSize) * 255);
         end
         
         % Get number of loops (the WHILE CPUTIME loop has more overhead):
         iTime = cputime;
         N = 0;
         while cputime - iTime < TestTime
            S = cell2mat(CStr);
            clear('S');
            N = N + 1;
         end
         N = max(2, round(N / (cputime - iTime) * SpeedTime));
         
         % The actual test:
         tic;
         for iN = 1:N
            S = [CStr{:}];
            clear('S');
         end
         NPerSec = N / (toc + eps);  % Loops per second
         PrintLoop(NPerSec);
         drawnow;  % Allow update of external events
      end
      fprintf('  loops/sec\n');
      
      % *********
      fprintf('    Cell2Vec(C)        ');
      for iCellSize = 1:length(CellSize)
         aCellSize = CellSize{iCellSize};
         if rejectLargeData && ...
               prod(aDataSize) >= 10000 && prod(aCellSize) >= 10000
            PrintLoop([]);
            continue;
         end
         
         CStr = cell(aCellSize);
         for iC = 1:prod(aCellSize)
            CStr{iC} = feval(aType, rand(aDataSize) * 255);
         end
         
         % Get number of loops (the WHILE CPUTIME loop has more overhead):
         iTime = cputime;
         N = 0;
         while cputime - iTime < TestTime
            S = Cell2Vec(CStr);
            clear('S');
            N = N + 1;
         end
         N = max(2, round(N / (cputime - iTime) * SpeedTime));
         
         % The actual test:
         tic;
         for iN = 1:N
            S = Cell2Vec(CStr);
            clear('S');
         end
         NPerSec = N / (toc + eps);  % Loops per second
         PrintLoop(NPerSec);
         drawnow;  % Allow update of external events
      end
      fprintf('  loops/sec\n');
   end  % for iDataSize
   fprintf('\n');
end  % for iType

disp('Cell2Vec passed the tests.');

% return;

% ******************************************************************************
function PrintLoop(N)

if isempty(N)
   fprintf('     slow');
elseif N > 10
   fprintf('  %7.0f', N);
else
   fprintf('  %7.1f', N);
end

% return;

