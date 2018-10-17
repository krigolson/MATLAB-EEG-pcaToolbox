on run argv
set a to item 1 of argv
set b to item 2 of argv
set c to item 3 of argv

set this_file to POSIX file a
         tell application "Finder"
            set file type of file this_file to b
            set creator type of file this_file to c
         end tell
end run