<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D tecplot output with a range restriction</description>
    <executable>FieldConvert</executable>
    <parameters> -r -1,1,-1,1 -e bfs_tg.xml bfs_tg.fld bfs_tg.dat</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="bfs_tg.dat">
                <sha1>7d8fb6d1df918e7218f8d5705742e67a69fd51bc</sha1>
             </file>
         </metric>
    </metrics>
</test>

