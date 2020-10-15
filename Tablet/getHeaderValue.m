function [value,lineNr]=getHeaderValue(headerData,variableName)
% getHeaderValue  Get value from text header of datafile by variable name

% Jason Gallivan
% Modified by Anouk de Brouwer

firstLetter = variableName(1);
nLines = size(headerData,1);
varNameLength = size(variableName,2);

lineNr = 1;
done = 0;
while done==0
    lineText = headerData{lineNr};
    % search line of text for 'VariableName' followed by a space
    if strncmpi(lineText,firstLetter,1); % first letter only to speed it up
        if strncmpi(lineText,variableName,varNameLength) && isspace(lineText(varNameLength+1))
            value = lineText(varNameLength+2:end);
            if all(isstrprop(value,'digit') | isstrprop(value,'punct')) % check if char is numeric
                numValue = str2double(value);                           % convert string to numeric
                if ~isnan(numValue)                                     % use numeric value only if it is not NaN
                    value = numValue;
                end
            end
            done = 1;
        end
    end
    
    % increase line number if we're not at the end of the file yet
    if lineNr >= nLines
        % variable name not found
        value = NaN;
        lineNr = NaN;
        done = 1;
    end
    lineNr = lineNr + 1;
end

lineNr = lineNr - 1;
