function [fid,fname,pathname]=get_fid(perm, mask, title)
%function [fid,fname,pathname]=get_fid(perm, mask, title)
%
%This function prompts for a filename and pathname using the ui.
%It check first to see if a non-null string is returned from the ui.
%If a valid filename is entered, the function attempts to open the
%file with the specified permission, and returns a fid and 
%pathname+filename string.
%
%Argument perm is a file permission string (c.f. fopen). Arguments
%mask and title are passed to uigetfile (c.f.).
%

% Modification history:
% 01/22/96 PJ Added pathname to argout list

if nargin < 1
    msg{1}='Function requires a file permision string as the first argument.';
    [msg]=ep_errorMsg(msg);
    return
elseif nargin == 1
    mask='*.*';
    title='Open File:';
elseif nargin == 2
    title='Open File:';
end

[fname, pathname] = uigetfile(mask,title);
if isempty(fname)
    fid= -1;
    msg{1}='No filename selected. You have to click on a name';
    [msg]=ep_errorMsg(msg);
    return
end

fname = [pathname fname];

[fid, message]=fopen(fname,perm);
if fid == -1
    disp(message)
end




