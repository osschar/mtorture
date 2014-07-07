#!/bin/bash
################################################################################
# Core helpers
################################################################################

function after_mic_reset()
{
    scp /opt/intel/composer_xe_2013_sp1.0.080/compiler/lib/mic/libiomp5.so root@mic0:/usr/lib64
    scp /opt/intel/composer_xe_2013_sp1.0.080/compiler/lib/mic/libiomp5.so root@mic1:/usr/lib64
}

### Turbo setting ###
#
# Turbo is (it seems):
# - bit 38 of 0x1a0: IA32_MISC_TURBO_EN  0x1a0 *** Inverse meaning, 0 means on ***
# - bit 32 of 0x199: IA32_PERF_TURBO_DIS 0x199
# The following doesn't really change anything ... there must be something
# else somewhere,

function echo_turbo()
{
    local t1
    local t2
    echo "Echo turbo (1-enabled, 0-disabled):"
    echo "Core | IA32_MISC_TURBO_EN | IA32_PERF_TURBO_DIS"
    for i in $(seq 0 23); do
        t1=`rdmsr -p $i 0x1a0 -f 38:38`
        t1=$(( 1 - $t1 ))
        t2=`rdmsr -p $i 0x199 -f 32:32`
        printf " %2d  | %1d                  | %1d\n" $i $t1 $t2
    done
}

function enable_turbo()
{
    local r

    for i in $(seq 0 23); do
        r=0x`rdmsr -p $i 0x1a0`
        r=`printf "0x%llx" $((~(1 << 38) & $r))`
        wrmsr -p $i 0x1a0 $r

        r=0x`rdmsr -p $i 0x199`
        r=`printf "0x%llx" $(( (1 << 32) | $r))`
        wrmsr -p $i 0x199 $r
    done
}

function disable_turbo()
{
    local r

    for i in $(seq 0 23); do
        r=0x`rdmsr -p $i 0x1a0`
        r=`printf "0x%llx" $(( (1 << 38) | $r))`
        wrmsr -p $i 0x1a0 $r

        r=0x`rdmsr -p $i 0x199`
        r=`printf "0x%llx" $((~(1 << 32) & $r))`
        wrmsr -p $i 0x199 $r
    done
}

################################################################################
# User / key management
################################################################################

function merge_auth_keys_for_root()
{
    rm -f authorized_keys
    for user in matevz; do
	cat /home/$user/.ssh/authorized_keys >> authorized_keys
    done

    scp authorized_keys mic0:.ssh/
    scp authorized_keys mic1:.ssh/
}

function clone_ssh_dir()
{
    local u=$1

    if [ -z "$u" ]; then
        echo "Tell me the user ..."
        return
    elif [ ! -e "/home/$u/" ]; then
        echo "User does not seem to exist ..."
        return
    fi

    for m in root@mic0 root@mic1; do
        ssh $m adduser -D $u
	scp -r /home/$u/.ssh $m:/home/$u
	ssh $m chown -R $u:$u /home/$u
	ssh $m passwd -u $u
    done
}

function clone_ssh_dirs()
{
    # for u in `ls /home`; do
    for u in matevz dsr; do
	clone_ssh_dir $u
    done
}


function delete_matevz()
{
    for m in root@mic0 root@mic1; do
	ssh $m deluser matevz
	ssh $m rm -rf /home/matevz
    done
}	    

function create_user()
{
    local u=$1

    if [ -z "$u" ]; then
	echo "Tell me the user ..."
	return
    elif [ -e "/home/$u/" ]; then
	echo "User already exists ..."
	return
    fi
    
    pushd /tmp

    useradd "$u"
    cd "/home/$u"
    chmod 0755 .
    mkdir .ssh
    touch .ssh/authorized_keys
    chown -R "$u:$u" .ssh
    chmod 0700 .ssh
    chmod 0600 .ssh/authorized_keys

    popd

    echo "cat > /home/$u/.ssh/authorized_keys"
    echo "clone_ssh_dir $u"
}
