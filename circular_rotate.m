function out = circular_rotate(A, maxIndex)
% [~, maxIndex] = max(A);
middleIndex = ceil(length(A) / 2);
shiftAmount = middleIndex - maxIndex;
A_rotated = circshift(A, shiftAmount);
out = [A_rotated(end), A_rotated];
end