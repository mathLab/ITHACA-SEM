<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D vorticity output </description>
    <executable>FieldConvert</executable>
    <parameters> -m vorticity bfs_tg.xml bfs_tg.fld bfs_tg_vort.fld</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="bfs_tg_vort.fld">
                <sha1>0ea257823268f80ad31b77f3ed89dfd2c818b39d</sha1> 
            </file>
         </metric>
    </metrics>
</test>

