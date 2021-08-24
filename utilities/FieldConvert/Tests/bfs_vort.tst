<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> Process 2D vorticity output </description>
    <executable>FieldConvert</executable>
    <parameters> -f -m vorticity -e bfs_tg.xml bfs_tg.fld bfs_tg_vort.fld</parameters>
    <files>
        <file description="Session File">bfs_tg.xml</file>
	<file description="Session File">bfs_tg.fld</file>
    </files>
     <metrics>
        <metric type="L2" id="1">
            <value variable="x" tolerance="1e-6">289.252</value>
            <value variable="y" tolerance="1e-6">6.0553</value>
            <value variable="u" tolerance="1e-6">4.6773</value>
            <value variable="v" tolerance="1e-4">0.172191</value>
            <value variable="p" tolerance="1e-6">0.359627</value>
            <value variable="W_z" tolerance="1e-6">10.8071</value>
        </metric>
    </metrics>
</test>

