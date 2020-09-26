<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> NekMesh with Spherigons and variable Boundary Layer </description>
    <executable>NekMesh</executable>
    <parameters> -m spherigon:surf=10:surf=13 -m spherigon:surf=8:surf=9 -m bl:surf=3,10,13:layers=4:r="1.7*( 1-x/0.3 )+1":nq=7 -m bl:surf=2,8,9:layers=4:r="1.7*(1-(x-0.27)/0.078)+1":nq=7 -m jac:list spherigon_bl_straight_rw.dat spherigon_bl_straight_rw-out.xml:xml:test </parameters>
    <files>
        <file description="Input File">spherigon_bl_straight_rw.dat</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
