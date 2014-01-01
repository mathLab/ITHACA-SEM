<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D vorticity output </description>
    <executable>FieldConvert</executable>
    <parameters> -r -1,1,-1,1 -m vorticity bfs_tg.xml bfs_tg.fld bfs_tg_vort.fld</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="bfs_tg_vort.fld">
                <sha1>2f8703f27a9d1b717caa97d796742d0d8eec8632</sha1> 
            </file>
         </metric>
    </metrics>
</test>

