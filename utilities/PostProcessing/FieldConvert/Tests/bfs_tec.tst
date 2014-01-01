<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D tecplot output </description>
    <executable>FieldConvert</executable>
    <parameters> bfs_tg.xml bfs_tg.fld bfs_tg.dat</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="bfs_tg.dat">
                <sha1>3a02738c68399ff375a8ba44ec1c8aa8ef5ea95d</sha1> 
            </file>
         </metric>
    </metrics>
</test>

