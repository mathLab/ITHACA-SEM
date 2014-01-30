<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D tecplot output </description>
    <executable>FieldConvert</executable>
    <parameters> -e bfs_tg.xml bfs_tg.fld bfs_tg.dat</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="file" id="1">
            <file filename="bfs_tg.dat">
                <sha1>2a8aca916f8db040c2a8c2ced64cc3298ddfac04</sha1> 
            </file>
         </metric>
    </metrics>
</test>

