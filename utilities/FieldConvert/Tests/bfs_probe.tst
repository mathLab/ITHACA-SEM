<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Probe line of 2D points from 2D file </description>
    <executable>FieldConvert</executable>
    <parameters> -f -e -m interppoints:fromxml=bfs_tg.xml:fromfld=bfs_tg.fld:topts=bfs_probe.pts bfs_probe.dat</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
	<file description="Session File">bfs_probe.pts</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">0.5</value>
            <value variable="y" tolerance="1e-6">0.5</value>
            <value variable="z" tolerance="1e-6">0.0</value>
            <value variable="u" tolerance="1e-6">0.577227</value>
            <value variable="v" tolerance="1e-6">0.00282065</value>
            <value variable="p" tolerance="1e-6">0.0333815</value>
        </metric>
    </metrics>
</test>

